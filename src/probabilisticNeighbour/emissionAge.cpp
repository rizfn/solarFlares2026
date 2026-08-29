// Tracks how long each spot has existed and how far it has travelled since injection,
// and reports both as a function of the emission size it eventually produces. Tests the
// picture that small emissions are spots injected straight into a wrong-sign domain
// (dying young and small) while large ones are spots injected into their own domain
// that grow while random-walking until they reach a wall.
// Splits the emission spectrum by the mass ratio of the annihilating pair, to test
// whether the two-slope structure at intermediate p is a mixture of two channels:
// unequal-mass annihilation (a small spot dying inside a domain) and comparable-mass
// annihilation (two domain-scale spots meeting). Also logs the branching ratio q
// resolved by mass ratio, which is the quantity the sub-leading balance is sensitive to.
// usage: ./emissionChannels L rho steps p seed [outDir]

#include <cmath>
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

constexpr long long DEFAULT_L = 128;
constexpr double DEFAULT_RHO = 0.2;
constexpr long long DEFAULT_STEPS = 100000;
constexpr double DEFAULT_P = 1.0;
constexpr int DEFAULT_NSNAP = 10;
constexpr int RECORD_INTERVAL = 100;    // sweeps between spot-histogram samples
constexpr int MAX_INTERVAL = 100;       // sweeps between largest-spot samples

constexpr double R_UNEQ = 0.10;   // "one partner is negligible"
constexpr double R_COMP = 0.25;   // "both partners carry comparable mass"
constexpr long long Q_MINMASS = 30;   // ignore the injection scale when tallying q
constexpr int NQ = 10;                // mass-ratio bins for q
long long qC[NQ] = {0}, qA[NQ] = {0};
long long curStep = 0;

constexpr int NAGE = 40;              // log bins in emission size, 6 per decade
double ageSum[NAGE] = {0}, dispSum[NAGE] = {0}, disp2Sum[NAGE] = {0};
long long ageN[NAGE] = {0};
std::vector<long long> bstep;         // sweep at which the spot now here was injected
std::vector<int> bx, by;              // and where
std::vector<char> bmode;              // 1 = injected next to its own sign, 0 = at random
std::vector<long long> ncoll;         // collisions this spot has already survived

inline int agebin(long long v)
{
    int b = (int)(6.0 * std::log10((double)v + 0.5));
    return b < 0 ? 0 : (b >= NAGE ? NAGE - 1 : b);
}

inline void born(int loc, int L, long long step, char mode)
{
    bstep[loc] = step; bx[loc] = loc % L; by[loc] = loc / L; bmode[loc] = mode; ncoll[loc] = 0;
}

inline void carry(int dst, int src_)
{
    bstep[dst] = bstep[src_]; bx[dst] = bx[src_]; by[dst] = by[src_]; bmode[dst] = bmode[src_];
    ncoll[dst] = ncoll[src_];
}

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

Hist *emisOwn = nullptr, *emisRand = nullptr;   // emissions split by that origin
constexpr int NLIFE = 40;             // log bins in lifetime, 6 per decade
long long lifeOwn[NLIFE] = {0}, lifeRand[NLIFE] = {0};
long long joint[NLIFE][NAGE] = {{0}};   // (lifetime, emitted size) at each annihilation
long long qAgeC[NLIFE] = {0}, qAgeA[NLIFE] = {0};
long long qNC[NLIFE] = {0}, qNA[NLIFE] = {0};      // merge vs annihilate, by that count
double ncSizeSum[NLIFE] = {0}; long long ncSizeN[NLIFE] = {0};
long long jointN[NLIFE][NAGE] = {{0}};   // (collisions survived, emitted size)  // merge vs annihilate, by the age of
                                                  // the smaller partner, over ALL collisions

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

// neighbour rule: place +1 next to a +, -1 next to a -
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
                if (s[nb] == 0) { s[nb] = 1; addFilled(nb, f, pos); born(nb, L, curStep, 1); posF = true; break; }
        }
        else if (s[loc] < 0 && !negF)
        {
            for (int nb : neighbours(loc, L))
                if (s[nb] == 0) { s[nb] = -1; addFilled(nb, f, pos); born(nb, L, curStep, 1); negF = true; break; }
        }
    }
    std::uniform_int_distribution<> dl(0, L * L - 1);
    if (!posF) { int r; do { r = dl(gen); } while (s[r] != 0); s[r] = 1; addFilled(r, f, pos); born(r, L, curStep, 0); }
    if (!negF) { int r; do { r = dl(gen); } while (s[r] != 0); s[r] = -1; addFilled(r, f, pos); born(r, L, curStep, 0); }
}

// random rule: + and - land on two independent empty sites
void addRandomPair(std::vector<long long> &s, int L, std::vector<int> &f, std::vector<int> &pos)
{
    std::uniform_int_distribution<> dl(0, L * L - 1);
    int a; do { a = dl(gen); } while (s[a] != 0);
    s[a] = 1; addFilled(a, f, pos); born(a, L, curStep, 0);
    int b; do { b = dl(gen); } while (s[b] != 0);
    s[b] = -1; addFilled(b, f, pos); born(b, L, curStep, 0);
}

void update(std::vector<long long> &s, int L, int N, double p, bool record, Hist &emis,
            Hist &emisU, Hist &emisC, Hist &coll, std::vector<int> &f, std::vector<int> &pos)
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
            double r = (double)lo / (double)hi;
            coll.add(lo);   // min mass over ALL collisions, the unweighted baseline
            {   // and the same collisions sorted by the age of the smaller partner
                int sm2 = (am <= ad) ? loc : dst;
                long long g = curStep - bstep[sm2];
                int gb = agebin(g > 0 ? g : 1);
                if (val * dv > 0) ++qAgeC[gb]; else ++qAgeA[gb];
                int nb2 = agebin(ncoll[sm2] + 1);
                if (val * dv > 0) ++qNC[nb2]; else ++qNA[nb2];
                if (val * dv < 0)
                {
                    ncSizeSum[nb2] += (double)lo; ++ncSizeN[nb2];
                    ++jointN[nb2][agebin(lo)];
                }
            }
            if (val * dv < 0)
            {
                emis.add(lo);
                if (r < R_UNEQ) emisU.add(lo);
                if (r > R_COMP) emisC.add(lo);
                // the emitted mass is the smaller partner: report ITS history
                int sm = (am <= ad) ? loc : dst;
                int cx = sm % L, cy = sm / L;
                int ddx = std::abs(cx - bx[sm]); ddx = std::min(ddx, L - ddx);
                int ddy = std::abs(cy - by[sm]); ddy = std::min(ddy, L - ddy);
                double d = std::sqrt((double)(ddx * ddx + ddy * ddy));
                int b2 = agebin(lo);
                ageSum[b2] += (double)(curStep - bstep[sm]);
                dispSum[b2] += d; disp2Sum[b2] += d * d; ++ageN[b2];
                long long life = curStep - bstep[sm];
                int lb = agebin(life > 0 ? life : 1);
                if (bmode[sm]) { emisOwn->add(lo); ++lifeOwn[lb]; }
                else           { emisRand->add(lo); ++lifeRand[lb]; }
                ++joint[lb][b2];
            }
            if (hi > Q_MINMASS)
            {
                int b = std::min(NQ - 1, (int)(r * NQ));
                if (val * dv > 0) ++qC[b]; else ++qA[b];
            }
        }
        if (val == -dv) removeLoc(dst, f, pos);
        // the survivor keeps the heavier partner's history
        if (std::llabs(val) > std::llabs(dv)) carry(dst, loc);
        ++ncoll[dst];
    }
    else { addFilled(dst, f, pos); carry(dst, loc); }
    s[dst] += val; s[loc] -= val;

    if ((int)f.size() < N)
    {
        if (dr(gen) < p) addNeighbourPair(s, L, f, pos);
        else addRandomPair(s, L, f, pos);
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

void run(int L, double rho, long long steps, double p,
         std::ofstream &spotF, std::ofstream &emisF, std::ofstream &emisUF,
         std::ofstream &emisCF, std::ofstream &collF, std::ofstream &qF, std::ofstream &ageF,
         std::ofstream &ownF, std::ofstream &randF, std::ofstream &lifeF, std::ofstream &jointF,
         std::ofstream &qAgeF, std::ofstream &ncF,
         std::ofstream &jnF)
{
    int N = (int)(((long long)L * L * rho) / 2) * 2;
    std::vector<long long> s((size_t)L * L, 0);
    std::fill(s.begin(), s.begin() + N / 2, 1);
    std::fill(s.begin() + N / 2, s.begin() + N, -1);
    std::shuffle(s.begin(), s.end(), gen);

    std::vector<int> f; f.reserve(2 * N + 8);
    std::vector<int> pos((size_t)L * L, -1);
    bstep.assign((size_t)L * L, 0); bx.assign((size_t)L * L, 0); by.assign((size_t)L * L, 0);
    bmode.assign((size_t)L * L, 0); ncoll.assign((size_t)L * L, 0);
    for (long long i = 0; i < (long long)L * L; ++i) if (s[i]) born((int)i, L, 0, 0);
    for (long long i = 0; i < (long long)L * L; ++i) if (s[i]) addFilled((int)i, f, pos);

    Hist spot, emis, emisU, emisC, coll, eOwn, eRand;
    emisOwn = &eOwn; emisRand = &eRand;
    long long rec = steps / 2;   // discard the first half as transient
    for (long long step = 0; step < steps; ++step)
    {
        bool record = step >= rec;
        curStep = step;
        for (long long i = 0; i < (long long)L * L; ++i) update(s, L, N, p, record, emis, emisU, emisC, coll, f, pos);
        if (record && step % RECORD_INTERVAL == 0)
            for (long long v : s) if (v) spot.add(std::llabs(v));
    }
    spot.write(spotF); emis.write(emisF); emisU.write(emisUF); emisC.write(emisCF); coll.write(collF);
    eOwn.write(ownF); eRand.write(randF);
    jnF << "# nCollisions\tsize\tcount\n";
    for (int i = 0; i < NLIFE; ++i)
        for (int j = 0; j < NAGE; ++j)
            if (jointN[i][j] > 0)
                jnF << std::pow(10.0, (i + 0.5) / 6.0) << "\t"
                    << std::pow(10.0, (j + 0.5) / 6.0) << "\t" << jointN[i][j] << "\n";
    ncF << "# nCollisions\tnMerge\tnAnnihilate\tmean_emitted_size\n";
    for (int b = 0; b < NLIFE; ++b)
        if (qNC[b] + qNA[b] > 200)
            ncF << std::pow(10.0, (b + 0.5) / 6.0) << "\t" << qNC[b] << "\t" << qNA[b]
                << "\t" << (ncSizeN[b] ? ncSizeSum[b] / ncSizeN[b] : 0.0) << "\n";
    qAgeF << "# age\tnMerge\tnAnnihilate\n";
    for (int b = 0; b < NLIFE; ++b)
        if (qAgeC[b] + qAgeA[b] > 200)
            qAgeF << std::pow(10.0, (b + 0.5) / 6.0) << "\t" << qAgeC[b] << "\t"
                  << qAgeA[b] << "\n";
    jointF << "# lifetime\tsize\tcount\n";
    for (int i = 0; i < NLIFE; ++i)
        for (int j = 0; j < NAGE; ++j)
            if (joint[i][j] > 0)
                jointF << std::pow(10.0, (i + 0.5) / 6.0) << "\t"
                       << std::pow(10.0, (j + 0.5) / 6.0) << "\t" << joint[i][j] << "\n";
    lifeF << "# lifetime\tn_own\tn_random\n";
    for (int b = 0; b < NLIFE; ++b)
        if (lifeOwn[b] + lifeRand[b] > 50)
            lifeF << std::pow(10.0, (b + 0.5) / 6.0) << "\t" << lifeOwn[b] << "\t"
                  << lifeRand[b] << "\n";
    ageF << "# s\tcount\tmean_age_sweeps\tmean_disp\trms_disp\n";
    for (int b = 0; b < NAGE; ++b)
        if (ageN[b] > 50)
            ageF << std::pow(10.0, (b + 0.5) / 6.0) << "\t" << ageN[b] << "\t"
                 << ageSum[b] / ageN[b] << "\t" << dispSum[b] / ageN[b] << "\t"
                 << std::sqrt(disp2Sum[b] / ageN[b]) << "\n";
    qF << "# ratio_lo\tnCoag\tnAnnih\tq\n";
    for (int b = 0; b < NQ; ++b)
    {
        double tot = (double)(qC[b] + qA[b]);
        qF << (double)b / NQ << "\t" << qC[b] << "\t" << qA[b] << "\t"
           << (tot > 0 ? qC[b] / tot : 0.0) << "\n";
    }
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

    // same filenames as the laptop version, so plots.py reads them unchanged
    std::ofstream spotF(outDir + "/chSpot_" + tag.str() + ".tsv");
    std::ofstream emisF(outDir + "/chEmisAll_" + tag.str() + ".tsv");
    std::ofstream emisUF(outDir + "/chEmisUneq_" + tag.str() + ".tsv");
    std::ofstream emisCF(outDir + "/chEmisComp_" + tag.str() + ".tsv");
    std::ofstream collF(outDir + "/chCollMin_" + tag.str() + ".tsv");
    std::ofstream ownF(outDir + "/chEmisOwn_" + tag.str() + ".tsv");
    std::ofstream randF(outDir + "/chEmisRand_" + tag.str() + ".tsv");
    std::ofstream jnF(outDir + "/chJointN_" + tag.str() + ".tsv");
    std::ofstream ncF(outDir + "/chNcoll_" + tag.str() + ".tsv");
    std::ofstream qAgeF(outDir + "/chQage_" + tag.str() + ".tsv");
    std::ofstream jointF(outDir + "/chJoint_" + tag.str() + ".tsv");
    std::ofstream lifeF(outDir + "/chLife_" + tag.str() + ".tsv");
    std::ofstream ageF(outDir + "/chAge_" + tag.str() + ".tsv");
    std::ofstream qF(outDir + "/chQratio_" + tag.str() + ".tsv");
    run((int)L, rho, steps, p, spotF, emisF, emisUF, emisCF, collF, qF, ageF, ownF, randF, lifeF, jointF, qAgeF, ncF, jnF);
    return 0;
}
