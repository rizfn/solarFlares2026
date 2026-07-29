// Instrumented copy of probabilisticNeighbour.cpp, used to test the assumption behind
// tau_s = 2 tau_m - 1 (that the two annihilating masses are independent draws from n(m)):
//  - dumps the (m1,m2) mass pair of every recorded annihilation event (subsampled)
//  - dumps a time series of max mass / total |mass| to check stationarity
//  - splits the spot histogram into early/late halves of the recording window
// usage: ./collisionPairs L rho steps p seed outDir [pairSkip]
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

constexpr int RECORD_INTERVAL = 100;

std::random_device rd;
std::mt19937 gen(rd());

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
    pos[loc] = (int)f.size(); f.push_back(loc);
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

long long pairSkip = 0, pairCtr = 0;

void update(std::vector<long long> &s, int L, int N, double p, bool record, Hist &emis,
            std::vector<int> &f, std::vector<int> &pos, std::ofstream &pairF)
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
        if (record && val * dv < 0)
        {
            emis.add(std::min(std::llabs(val), std::llabs(dv)));
            // mover mass first, target mass second: the two are distinguishable
            if (pairCtr++ % pairSkip == 0)
                pairF << std::llabs(val) << "\t" << std::llabs(dv) << "\n";
        }
        if (val == -dv) removeLoc(dst, f, pos);
    }
    else addFilled(dst, f, pos);
    s[dst] += val; s[loc] -= val;

    if ((int)f.size() < N)
    {
        if (dr(gen) < p) addNeighbourPair(s, L, f, pos);
        else addBipole(s, L, f, pos);
    }
}

int main(int argc, char *argv[])
{
    long long L = std::stoll(argv[1]), steps = std::stoll(argv[3]);
    double rho = std::stod(argv[2]), p = std::stod(argv[4]);
    unsigned seed = (unsigned)std::stoul(argv[5]);
    std::string outDir = argv[6];
    pairSkip = argc > 7 ? std::stoll(argv[7]) : 200;
    gen.seed(seed);
    std::filesystem::create_directories(outDir);
    std::ostringstream tag;
    tag << "L_" << L << "_rho_" << rho << "_p_" << p << "_seed_" << seed;
    std::ofstream pairF(outDir + "/pairs_" + tag.str() + ".tsv");
    std::ofstream tsF(outDir + "/timeseries_" + tag.str() + ".tsv");
    std::ofstream hEarly(outDir + "/spotEarly_" + tag.str() + ".tsv");
    std::ofstream hLate(outDir + "/spotLate_" + tag.str() + ".tsv");
    std::ofstream emisF(outDir + "/emisDiag_" + tag.str() + ".tsv");
    pairF << "# m_mover\tm_target\n";
    tsF << "# step\tnspots\ttotal_abs_mass\tmax_mass\n";

    int N = (int)((L * L * rho) / 2) * 2;
    std::vector<long long> s(L * L, 0);
    std::fill(s.begin(), s.begin() + N / 2, 1);
    std::fill(s.begin() + N / 2, s.begin() + N, -1);
    std::shuffle(s.begin(), s.end(), gen);
    std::vector<int> f; f.reserve(2 * N + 8);
    std::vector<int> pos(L * L, -1);
    for (int i = 0; i < L * L; ++i) if (s[i]) addFilled(i, f, pos);

    Hist emis, early, late;
    long long rec = steps / 2, mid = rec + (steps - rec) / 2;
    for (long long step = 0; step < steps; ++step)
    {
        bool record = step >= rec;
        for (int i = 0; i < L * L; ++i) update(s, L, N, p, record, emis, f, pos, pairF);
        if (step % RECORD_INTERVAL == 0)
        {
            long long tot = 0, mx = 0;
            for (long long v : s) { tot += std::llabs(v); mx = std::max(mx, std::llabs(v)); }
            tsF << step << "\t" << f.size() << "\t" << tot << "\t" << mx << "\n";
            if (record)
                for (long long v : s) if (v) (step < mid ? early : late).add(std::llabs(v));
        }
    }
    early.write(hEarly); late.write(hLate); emis.write(emisF);
    return 0;
}
