// Server build of probabilisticNeighbour. Physics identical to
// src/probabilisticNeighbour/probabilisticNeighbour.cpp; only the output path,
// snapshot count and sampling intervals differ. See readme.md.
// usage: ./probabilisticNeighbourServer L rho steps p seed [nSnapshots] [outDir]

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

constexpr const char *DEFAULT_OUTDIR = "/nbi/home/rpw391/cell-disk/solarFlares/outputs";

constexpr long long DEFAULT_L = 128;
constexpr double DEFAULT_RHO = 0.2;
constexpr long long DEFAULT_STEPS = 100000;
constexpr double DEFAULT_P = 1.0;
constexpr int DEFAULT_NSNAP = 10;
constexpr int RECORD_INTERVAL = 100;    // sweeps between spot-histogram samples
constexpr int MAX_INTERVAL = 100;       // sweeps between largest-spot samples

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

// random rule: + and - land on two independent empty sites
void addRandomPair(std::vector<long long> &s, int L, std::vector<int> &f, std::vector<int> &pos)
{
    std::uniform_int_distribution<> dl(0, L * L - 1);
    int a; do { a = dl(gen); } while (s[a] != 0);
    s[a] = 1; addFilled(a, f, pos);
    int b; do { b = dl(gen); } while (s[b] != 0);
    s[b] = -1; addFilled(b, f, pos);
}

void update(std::vector<long long> &s, int L, int N, double p, bool record, Hist &emis,
            std::vector<int> &f, std::vector<int> &pos)
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
        if (record && val * dv < 0) emis.add(std::min(std::llabs(val), std::llabs(dv)));
        if (val == -dv) removeLoc(dst, f, pos);
    }
    else addFilled(dst, f, pos);
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

void run(int L, double rho, long long steps, double p, int nSnap,
         std::ofstream &spotF, std::ofstream &emisF, std::ofstream &snapF,
         std::ofstream &maxF)
{
    int N = (int)(((long long)L * L * rho) / 2) * 2;
    std::vector<long long> s((size_t)L * L, 0);
    std::fill(s.begin(), s.begin() + N / 2, 1);
    std::fill(s.begin() + N / 2, s.begin() + N, -1);
    std::shuffle(s.begin(), s.end(), gen);

    std::vector<int> f; f.reserve(2 * N + 8);
    std::vector<int> pos((size_t)L * L, -1);
    for (long long i = 0; i < (long long)L * L; ++i) if (s[i]) addFilled((int)i, f, pos);

    Hist spot, emis;
    long long rec = steps / 2;   // discard the first half as transient
    long long snapEvery = std::max<long long>(1, (steps - rec) / std::max(nSnap, 1));
    maxF << "# step\tmaxAbs\tsumAbs\tnSpots\n";
    for (long long step = 0; step < steps; ++step)
    {
        bool record = step >= rec;
        for (long long i = 0; i < (long long)L * L; ++i) update(s, L, N, p, record, emis, f, pos);
        if (step % MAX_INTERVAL == 0)
        {
            long long mx = 0, tot = 0;
            for (int loc : f) { long long a = std::llabs(s[loc]); if (a > mx) mx = a; tot += a; }
            maxF << step << "\t" << mx << "\t" << tot << "\t" << f.size() << "\n";
        }
        if (record && step % RECORD_INTERVAL == 0)
            for (long long v : s) if (v) spot.add(std::llabs(v));
        if (record && nSnap > 0 && (step - rec) % snapEvery == 0)
            writeSnapshot(snapF, s, L);
    }
    spot.write(spotF); emis.write(emisF);
}

int main(int argc, char *argv[])
{
    long long L = DEFAULT_L, steps = DEFAULT_STEPS;
    double rho = DEFAULT_RHO, p = DEFAULT_P;
    unsigned seed = 1;
    int nSnap = DEFAULT_NSNAP;
    std::string outDir = DEFAULT_OUTDIR;
    if (argc > 1) L = std::stoll(argv[1]);
    if (argc > 2) rho = std::stod(argv[2]);
    if (argc > 3) steps = std::stoll(argv[3]);
    if (argc > 4) p = std::stod(argv[4]);
    if (argc > 5) seed = (unsigned)std::stoul(argv[5]);
    if (argc > 6) nSnap = std::stoi(argv[6]);
    if (argc > 7) outDir = argv[7];
    gen.seed(seed);

    std::filesystem::create_directories(outDir);
    std::ostringstream tag;
    tag << "L_" << L << "_rho_" << rho << "_p_" << p << "_seed_" << seed;

    // same filenames as the laptop version, so plots.py reads them unchanged
    std::ofstream spotF(outDir + "/spotSize_" + tag.str() + ".tsv");
    std::ofstream emisF(outDir + "/emission_" + tag.str() + ".tsv");
    std::ofstream snapF(outDir + "/snapshots_" + tag.str() + ".tsv");
    std::ofstream maxF(outDir + "/maxSpot_" + tag.str() + ".tsv");
    run((int)L, rho, steps, p, nSnap, spotF, emisF, snapF, maxF);
    return 0;
}
