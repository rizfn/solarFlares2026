// Injection correlation instead of injection sign-sorting.
//
// Same cascade as probabilisticNeighbour: spots random-walk, same signs coagulate
// (i,j -> i+j), opposite signs partially annihilate and emit s = min(i,j). The number
// of spots is held at N by injecting a +1/-1 pair whenever one is lost.
//
// The parameter is where that pair goes:
//   with probability r      the + and - are placed as lattice NEIGHBOURS (a bipole),
//   with probability (1-r)  they are placed at two INDEPENDENT random empty sites.
//
// r = 1 injects a maximally anti-correlated pair, which tends to annihilate before it
// can find anything else, so the same-sign collision fraction q drops below 1/2 — the
// regime where the mean-field rate equation has no scale-free solution at all. r = 0
// injects with no local correlation. Whether 2D drift restores a power law where mean
// field forbids one is the question this sweep is for, so q is measured, not imposed.
//
// usage: ./addedBipole L rho steps r seed outDir
#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <sstream>
#include <filesystem>

constexpr long long DEFAULT_L = 128;
constexpr double DEFAULT_RHO = 0.2;
constexpr long long DEFAULT_STEPS = 100000;
constexpr double DEFAULT_R = 1.0;
constexpr int RECORD_INTERVAL = 100;
constexpr int NUM_SNAPSHOTS = 40;

std::random_device rd;
std::mt19937 gen(rd());

// Sparse power-law histogram: dense small bins, hash map for the rare large tail.
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

// bipole rule: the + and the - land on adjacent sites
void addBipole(std::vector<long long> &s, int L, std::vector<int> &f, std::vector<int> &pos)
{
    std::uniform_int_distribution<> dl(0, L * L - 1);
    int p; do { p = dl(gen); } while (s[p] != 0);
    auto nb = neighbours(p, L);
    std::shuffle(nb.begin(), nb.end(), gen);
    int m = -1;
    for (int c : nb) if (s[c] == 0) { m = c; break; }
    if (m < 0) do { m = dl(gen); } while (s[m] != 0);
    s[p] = 1; addFilled(p, f, pos); s[m] = -1; addFilled(m, f, pos);
}

// random rule: the + and the - land on two independent empty sites. The + is written
// before the - is drawn, so the two cannot collide.
void addRandomPair(std::vector<long long> &s, int L, std::vector<int> &f, std::vector<int> &pos)
{
    std::uniform_int_distribution<> dl(0, L * L - 1);
    int a; do { a = dl(gen); } while (s[a] != 0);
    s[a] = 1; addFilled(a, f, pos);
    int b; do { b = dl(gen); } while (s[b] != 0);
    s[b] = -1; addFilled(b, f, pos);
}

void update(std::vector<long long> &s, int L, int N, double r, bool record, Hist &emis,
            std::vector<int> &f, std::vector<int> &pos, long long &nCoag, long long &nAnnih)
{
    std::uniform_int_distribution<> di(0, (int)f.size() - 1);
    int idx = di(gen), loc = f[idx];
    long long val = s[loc];
    std::uniform_real_distribution<> dr(0.0, 1.0);
    int x = loc % L, y = loc / L, dst;
    double u = dr(gen);
    if (u < 0.25) dst = ((y - 1 + L) % L) * L + x;
    else if (u < 0.5) dst = ((y + 1) % L) * L + x;
    else if (u < 0.75) dst = y * L + (x - 1 + L) % L;
    else dst = y * L + (x + 1) % L;

    long long dv = s[dst];
    removeAt(idx, f, pos);
    if (dv != 0)
    {
        // branching ratio q of the kinetic theory, measured rather than imposed
        if (record) { if (val * dv > 0) ++nCoag; else ++nAnnih; }
        if (record && val * dv < 0) emis.add(std::min(std::llabs(val), std::llabs(dv)));
        if (val == -dv) removeLoc(dst, f, pos);
    }
    else addFilled(dst, f, pos);
    s[dst] += val; s[loc] -= val;

    if ((int)f.size() < N)
    {
        if (dr(gen) < r) addBipole(s, L, f, pos);
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

void run(int L, double rho, long long steps, double r,
         std::ofstream &spotF, std::ofstream &emisF, std::ofstream &snapF, std::ofstream &statF)
{
    int N = (int)((L * L * rho) / 2) * 2;
    std::vector<long long> s(L * L, 0);
    std::fill(s.begin(), s.begin() + N / 2, 1);
    std::fill(s.begin() + N / 2, s.begin() + N, -1);
    std::shuffle(s.begin(), s.end(), gen);

    std::vector<int> f; f.reserve(2 * N + 8);
    std::vector<int> pos(L * L, -1);
    for (int i = 0; i < L * L; ++i) if (s[i]) addFilled(i, f, pos);

    Hist spot, emis;
    long long nCoag = 0, nAnnih = 0;
    long long rec = steps / 2;
    long long snapEvery = std::max<long long>(1, (steps - rec) / NUM_SNAPSHOTS);
    for (long long step = 0; step < steps; ++step)
    {
        bool record = step >= rec;
        for (int i = 0; i < L * L; ++i) update(s, L, N, r, record, emis, f, pos, nCoag, nAnnih);
        if (record && step % RECORD_INTERVAL == 0)
            for (long long v : s) if (v) spot.add(std::llabs(v));
        if (record && (step - rec) % snapEvery == 0)
            writeSnapshot(snapF, s, L);
    }
    spot.write(spotF); emis.write(emisF);
    statF << "# coagulations\tannihilations\tq\n";
    statF << nCoag << "\t" << nAnnih << "\t"
          << (double)nCoag / (double)(nCoag + nAnnih) << "\n";
}

int main(int argc, char *argv[])
{
    long long L = DEFAULT_L, steps = DEFAULT_STEPS;
    double rho = DEFAULT_RHO, r = DEFAULT_R;
    unsigned seed = std::random_device{}();
    std::string outDir = "outputs";
    if (argc > 1) L = std::stoll(argv[1]);
    if (argc > 2) rho = std::stod(argv[2]);
    if (argc > 3) steps = std::stoll(argv[3]);
    if (argc > 4) r = std::stod(argv[4]);
    if (argc > 5) seed = (unsigned)std::stoul(argv[5]);
    if (argc > 6) outDir = argv[6];
    gen.seed(seed);

    std::filesystem::create_directories(outDir);
    std::ostringstream tag;
    tag << "L_" << L << "_rho_" << rho << "_r_" << r << "_seed_" << seed;
    std::ofstream spotF(outDir + "/spotSize_" + tag.str() + ".tsv");
    std::ofstream emisF(outDir + "/emission_" + tag.str() + ".tsv");
    std::ofstream snapF(outDir + "/snapshots_" + tag.str() + ".tsv");
    std::ofstream statF(outDir + "/stats_" + tag.str() + ".tsv");
    run(L, rho, steps, r, spotF, emisF, snapF, statF);
    return 0;
}
