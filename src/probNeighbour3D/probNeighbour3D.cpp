// 3D version of probabilisticNeighbour: spots random-walk on a cubic lattice, same
// signs coagulate (i,j -> i+j), opposite signs partially annihilate and emit min(i,j).
// The spot count is held fixed by injecting a +1/-1 pair whenever one is lost:
//   with probability p     the + is placed next to an existing +, the - next to a -
//   with probability (1-p) the + and - land on two independent random empty sites
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

constexpr long long DEFAULT_L = 32;
constexpr double DEFAULT_RHO = 0.2;
constexpr long long DEFAULT_STEPS = 50000;
constexpr double DEFAULT_P = 1.0;
constexpr int RECORD_INTERVAL = 100;
constexpr int NUM_SNAPSHOTS = 20;
constexpr int MAX_INTERVAL = 10;   // largest-spot sampling interval, in sweeps

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

// site index is x + L*y + L*L*z; six neighbours, periodic in all directions
std::array<int, 6> neighbours(int i, int L)
{
    int x = i % L, y = (i / L) % L, z = i / (L * L);
    int xy = L * L;
    return {z * xy + y * L + (x + 1) % L,
            z * xy + y * L + (x - 1 + L) % L,
            z * xy + ((y + 1) % L) * L + x,
            z * xy + ((y - 1 + L) % L) * L + x,
            ((z + 1) % L) * xy + y * L + x,
            ((z - 1 + L) % L) * xy + y * L + x};
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
    std::uniform_int_distribution<> dl(0, L * L * L - 1);
    if (!posF) { int r; do { r = dl(gen); } while (s[r] != 0); s[r] = 1; addFilled(r, f, pos); }
    if (!negF) { int r; do { r = dl(gen); } while (s[r] != 0); s[r] = -1; addFilled(r, f, pos); }
}

// random rule: the + and the - land on two independent empty sites. The + is written
// before the - is drawn, so the two cannot collide.
void addRandomPair(std::vector<long long> &s, int L, std::vector<int> &f, std::vector<int> &pos)
{
    std::uniform_int_distribution<> dl(0, L * L * L - 1);
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
    std::uniform_int_distribution<> dd(0, 5);
    int dst = neighbours(loc, L)[dd(gen)];

    long long dv = s[dst];
    removeAt(idx, f, pos);
    if (dv != 0)
    {
        if (record && val * dv < 0) emis.add(std::min(std::llabs(val), std::llabs(dv)));
        if (val == -dv) removeLoc(dst, f, pos);
    }
    else addFilled(dst, f, pos);
    s[dst] += val; s[loc] -= val;

    std::uniform_real_distribution<> dr(0.0, 1.0);
    if ((int)f.size() < N)
    {
        if (dr(gen) < p) addNeighbourPair(s, L, f, pos);
        else addRandomPair(s, L, f, pos);
    }
}

// dump the sign lattice (-1/0/1) as L*L rows of L values, z-major, blank line between
// snapshots. Python reshapes this into (L, L, L) indexed [z][y][x].
void writeSnapshot(std::ofstream &f, const std::vector<long long> &s, int L)
{
    for (int z = 0; z < L; ++z)
        for (int y = 0; y < L; ++y)
        {
            for (int x = 0; x < L; ++x)
            {
                long long v = s[z * L * L + y * L + x];
                f << (v > 0 ? 1 : (v < 0 ? -1 : 0)) << (x + 1 < L ? "\t" : "");
            }
            f << "\n";
        }
    f << "\n";
}

void run(int L, double rho, long long steps, double p,
         std::ofstream &spotF, std::ofstream &emisF, std::ofstream &snapF,
         std::ofstream &maxF)
{
    long long V = (long long)L * L * L;
    int N = (int)((V * rho) / 2) * 2;
    std::vector<long long> s(V, 0);
    std::fill(s.begin(), s.begin() + N / 2, 1);
    std::fill(s.begin() + N / 2, s.begin() + N, -1);
    std::shuffle(s.begin(), s.end(), gen);

    std::vector<int> f; f.reserve(2 * N + 8);
    std::vector<int> pos(V, -1);
    for (long long i = 0; i < V; ++i) if (s[i]) addFilled((int)i, f, pos);

    Hist spot, emis;
    long long rec = steps / 2;
    long long snapEvery = std::max<long long>(1, (steps - rec) / NUM_SNAPSHOTS);
    maxF << "# step\tmaxAbs\tsumAbs\tnSpots\n";
    for (long long step = 0; step < steps; ++step)
    {
        bool record = step >= rec;
        for (long long i = 0; i < V; ++i) update(s, L, N, p, record, emis, f, pos);
        if (step % MAX_INTERVAL == 0)
        {
            long long mx = 0, tot = 0;
            for (int loc : f) { long long a = std::llabs(s[loc]); if (a > mx) mx = a; tot += a; }
            maxF << step << "\t" << mx << "\t" << tot << "\t" << f.size() << "\n";
        }
        if (record && step % RECORD_INTERVAL == 0)
            for (long long v : s) if (v) spot.add(std::llabs(v));
        if (record && (step - rec) % snapEvery == 0)
            writeSnapshot(snapF, s, L);
    }
    spot.write(spotF); emis.write(emisF);
}

int main(int argc, char *argv[])
{
    long long L = DEFAULT_L, steps = DEFAULT_STEPS;
    double rho = DEFAULT_RHO, p = DEFAULT_P;
    unsigned seed = std::random_device{}();
    std::string outDir = "outputs";
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
    std::ofstream spotF(outDir + "/spotSize_" + tag.str() + ".tsv");
    std::ofstream emisF(outDir + "/emission_" + tag.str() + ".tsv");
    std::ofstream snapF(outDir + "/snapshots_" + tag.str() + ".tsv");
    std::ofstream maxF(outDir + "/maxSpot_" + tag.str() + ".tsv");
    run((int)L, rho, steps, p, spotF, emisF, snapF, maxF);
    return 0;
}
