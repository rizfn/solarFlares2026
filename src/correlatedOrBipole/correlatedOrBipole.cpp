// probabilisticNeighbour with the uncorrelated injection channel replaced by a bipole.
//
// The cascade is unchanged: spots random-walk on a periodic square lattice, same signs
// coagulate (i,j -> i+j), opposite signs partially annihilate and emit s = min(i,j), and
// the spot count is held fixed by injecting a +1/-1 pair whenever one is lost. Only the
// (1-p) channel differs:
//
//   probabilisticNeighbour     with prob (1-p): + and - land on two INDEPENDENT random
//                              empty sites, so each is dropped alone into whatever
//                              domain happens to be there
//   this model                 with prob (1-p): + and - land ADJACENT to each other at a
//                              random empty location -- a bipole, charge-neutral locally
//
// The purpose is to test whether the small-s branch of the emission spectrum is an
// injection artefact. In probabilisticNeighbour the annihilating fraction 1 - q(s) is
// U-shaped, and the falling branch is attributed to minority monomers dropped into
// wrong-sign domains, which die immediately and emit s ~ 1. A bipole cannot do that: the
// pair is its own nearest opposite sign, so it either self-annihilates at once (emitting
// s = 1 and restoring the count, a null event) or the two halves separate and join the
// bulk. If the attribution is right, the U should flatten and the crossover scale s*
// should weaken or vanish, while the exponent still slides with p.
//
// usage: ./correlatedOrBipole L rho steps p seed [outDir]

#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <string>

#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{ std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr); std::cout.tie(nullptr); return 0; }();

constexpr const char *DEFAULT_OUTDIR = "outputs";

constexpr long long DEFAULT_L = 512;
constexpr double DEFAULT_RHO = 0.2;
constexpr long long DEFAULT_STEPS = 100000;
constexpr double DEFAULT_P = 0.8;
constexpr int RECORD_INTERVAL = 100;   // sweeps between spot-histogram samples

constexpr long long Q_MINMASS = 30;    // ignore the injection scale when tallying q
constexpr int NQ = 10;                 // mass-ratio bins for q
long long qC[NQ] = {0}, qA[NQ] = {0};

long long nBipole = 0, nBipoleFail = 0;   // how often the bipole rule had to fall back

std::random_device rd;
std::mt19937 gen(rd());

// sparse power-law histogram: dense small bins, hash map for the rare large tail
struct Hist
{
    static constexpr long long LIMIT = 1LL << 22;
    std::vector<long long> small;
    std::unordered_map<long long, long long> large;
    Hist() : small(1024, 0) {}
    void add(long long v)
    {
        if (v >= LIMIT) { ++large[v]; return; }
        if (v >= (long long)small.size())
        {
            long long ns = small.size();
            while (ns <= v) ns *= 2;
            small.resize(std::min(ns, LIMIT), 0);
        }
        ++small[v];
    }
    void write(std::ofstream &f) const
    {
        f << "# size\tcount\n";
        for (long long i = 0; i < (long long)small.size(); ++i)
            if (small[i]) f << i << "\t" << small[i] << "\n";
        std::vector<long long> keys;
        for (const auto &kv : large) keys.push_back(kv.first);
        std::sort(keys.begin(), keys.end());
        for (long long k : keys) f << k << "\t" << large.at(k) << "\n";
    }
};

std::array<int, 4> neighbours(int i, int L)
{
    int x = i % L, y = i / L;
    return {((y - 1 + L) % L) * L + x, ((y + 1) % L) * L + x,
            y * L + (x - 1 + L) % L, y * L + (x + 1) % L};
}

void addFilled(int loc, std::vector<int> &f, std::vector<int> &pos)
{
    pos[loc] = (int)f.size();
    f.push_back(loc);
}
void removeAt(int idx, std::vector<int> &f, std::vector<int> &pos)
{
    int loc = f[idx], last = f.back();
    f[idx] = last; pos[last] = idx; f.pop_back(); pos[loc] = -1;
}
void removeLoc(int loc, std::vector<int> &f, std::vector<int> &pos)
{
    removeAt(pos[loc], f, pos);
}

// correlated rule, unchanged: place +1 next to a +, -1 next to a -
void addNeighbourPair(std::vector<long long> &s, int L, std::vector<int> &f, std::vector<int> &pos)
{
    bool posF = false, negF = false;
    int n = (int)f.size();
    for (int i = 0; i < n && !(posF && negF); ++i)
    {
        std::uniform_int_distribution<> d(i, n - 1);
        int j = d(gen), a = f[i], b = f[j];
        f[i] = b; f[j] = a; pos[b] = i; pos[a] = j;
        int loc = f[i];
        if (s[loc] > 0 && !posF)
        {
            for (int nb : neighbours(loc, L))
                if (s[nb] == 0) { s[nb] = 1; addFilled(nb, f, pos); posF = true; break; }
        }
        else if (s[loc] < 0 && !negF)
        {
            for (int nb : neighbours(loc, L))
                if (s[nb] == 0) { s[nb] = -1; addFilled(nb, f, pos); negF = true; break; }
        }
    }
    std::uniform_int_distribution<> dl(0, L * L - 1);
    if (!posF) { int r; do { r = dl(gen); } while (s[r] != 0); s[r] = 1; addFilled(r, f, pos); }
    if (!negF) { int r; do { r = dl(gen); } while (s[r] != 0); s[r] = -1; addFilled(r, f, pos); }
}

// bipole rule: + and - land on two ADJACENT empty sites, orientation random.
// Falls back to two independent sites only if no empty site has an empty neighbour,
// which at these densities essentially never happens; the count is reported so that
// "essentially never" is a measured statement rather than an assumption.
void addBipolePair(std::vector<long long> &s, int L, std::vector<int> &f, std::vector<int> &pos)
{
    std::uniform_int_distribution<> dl(0, L * L - 1);
    for (int attempt = 0; attempt < 64; ++attempt)
    {
        int a; do { a = dl(gen); } while (s[a] != 0);
        auto nb = neighbours(a, L);
        int cand[4], nc = 0;
        for (int k = 0; k < 4; ++k) if (s[nb[k]] == 0) cand[nc++] = nb[k];
        if (nc == 0) continue;
        std::uniform_int_distribution<> dc(0, nc - 1);
        int b = cand[dc(gen)];
        s[a] = 1;  addFilled(a, f, pos);
        s[b] = -1; addFilled(b, f, pos);
        ++nBipole;
        return;
    }
    ++nBipoleFail;
    int a; do { a = dl(gen); } while (s[a] != 0);
    s[a] = 1; addFilled(a, f, pos);
    int b; do { b = dl(gen); } while (s[b] != 0);
    s[b] = -1; addFilled(b, f, pos);
}

void update(std::vector<long long> &s, int L, int N, double p, bool record,
            Hist &emis, Hist &coll, std::vector<int> &f, std::vector<int> &pos)
{
    std::uniform_int_distribution<> di(0, (int)f.size() - 1);
    int idx = di(gen), loc = f[idx];
    long long val = s[loc];
    std::uniform_real_distribution<> dr(0.0, 1.0);
    int x = loc % L, y = loc / L, dst;
    double r = dr(gen);
    if (r < 0.25) dst = ((y - 1 + L) % L) * L + x;
    else if (r < 0.5) dst = ((y + 1) % L) * L + x;
    else if (r < 0.75) dst = y * L + (x - 1 + L) % L;
    else dst = y * L + (x + 1) % L;

    long long dv = s[dst];
    removeAt(idx, f, pos);
    if (dv != 0)
    {
        if (record)
        {
            long long am = std::llabs(val), ad = std::llabs(dv);
            long long lo = std::min(am, ad), hi = std::max(am, ad);
            coll.add(lo);                       // min mass over ALL collisions
            if (val * dv < 0) emis.add(lo);     // the annihilating subset
            if (hi > Q_MINMASS)
            {
                int b = std::min(NQ - 1, (int)((double)lo / (double)hi * NQ));
                if (val * dv > 0) ++qC[b]; else ++qA[b];
            }
        }
        if (val == -dv) removeLoc(dst, f, pos);
    }
    else addFilled(dst, f, pos);
    s[dst] += val; s[loc] -= val;

    if ((int)f.size() < N)
    {
        if (dr(gen) < p) addNeighbourPair(s, L, f, pos);
        else addBipolePair(s, L, f, pos);
    }
}

// dump the sign lattice (-1/0/1), snapshots separated by a blank line
void writeSnapshot(std::ofstream &f, const std::vector<long long> &s, int L)
{
    for (int y = 0; y < L; ++y)
    {
        for (int x = 0; x < L; ++x)
        {
            long long v = s[y * L + x];
            f << (v > 0 ? 1 : (v < 0 ? -1 : 0)) << (x + 1 < L ? "\t" : "");
        }
        f << "\n";
    }
    f << "\n";
}

void run(int L, double rho, long long steps, double p, std::ofstream &spotF,
         std::ofstream &emisF, std::ofstream &collF, std::ofstream &qF,
         std::ofstream &snapF)
{
    int N = (int)(((long long)L * L * rho) / 2) * 2;
    std::vector<long long> s((size_t)L * L, 0);
    std::fill(s.begin(), s.begin() + N / 2, 1);
    std::fill(s.begin() + N / 2, s.begin() + N, -1);
    std::shuffle(s.begin(), s.end(), gen);

    std::vector<int> f; f.reserve(2 * N + 8);
    std::vector<int> pos((size_t)L * L, -1);
    for (long long i = 0; i < (long long)L * L; ++i) if (s[i]) addFilled((int)i, f, pos);

    Hist spot, emis, coll;
    long long rec = steps / 2;   // discard the first half as transient
    for (long long step = 0; step < steps; ++step)
    {
        bool record = step >= rec;
        for (long long i = 0; i < (long long)L * L; ++i)
            update(s, L, N, p, record, emis, coll, f, pos);
        if (record && step % RECORD_INTERVAL == 0)
            for (long long v : s) if (v) spot.add(std::llabs(v));
    }
    writeSnapshot(snapF, s, L);
    spot.write(spotF); emis.write(emisF); coll.write(collF);

    qF << "# ratio_lo\tnCoag\tnAnnih\tq\n";
    for (int b = 0; b < NQ; ++b)
    {
        double tot = (double)(qC[b] + qA[b]);
        qF << (double)b / NQ << "\t" << qC[b] << "\t" << qA[b] << "\t"
           << (tot > 0 ? qC[b] / tot : 0.0) << "\n";
    }
    std::cerr << "bipole injections: " << nBipole << ", fallbacks: " << nBipoleFail << "\n";
}

int main(int argc, char *argv[])
{
    long long L = DEFAULT_L, steps = DEFAULT_STEPS;
    double rho = DEFAULT_RHO, p = DEFAULT_P;
    unsigned seed = 1;
    std::string outDir = DEFAULT_OUTDIR;
    if (argc > 1) L = std::stoll(argv[1]);
    if (argc > 2) rho = std::stod(argv[2]);
    if (argc > 3) steps = std::stoll(argv[3]);
    if (argc > 4) p = std::stod(argv[4]);
    if (argc > 5) seed = (unsigned)std::stoul(argv[5]);
    if (argc > 6) outDir = argv[6];
    gen.seed(seed);

    std::filesystem::create_directories(outDir);
    std::ostringstream tag;
    tag << "L_" << L << "_rho_" << rho << "_p_" << p << "_seed_" << seed;

    std::ofstream spotF(outDir + "/cbSpot_" + tag.str() + ".tsv");
    std::ofstream emisF(outDir + "/cbEmisAll_" + tag.str() + ".tsv");
    std::ofstream collF(outDir + "/cbCollMin_" + tag.str() + ".tsv");
    std::ofstream qF(outDir + "/cbQratio_" + tag.str() + ".tsv");
    std::ofstream snapF(outDir + "/cbSnapshot_" + tag.str() + ".tsv");
    run((int)L, rho, steps, p, spotF, emisF, collF, qF, snapF);
    return 0;
}
