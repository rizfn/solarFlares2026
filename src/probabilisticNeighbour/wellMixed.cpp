// Well-mixed signed coagulation-annihilation cascade: the Takayasu model with two
// signs and biased partner selection, on a complete graph instead of a lattice.
//
// The branching ratio q is imposed directly as a rule rather than emerging from
// spatial correlations, so this is the exact stochastic process whose deterministic
// limit is the rate equation integrated by meanField.py. That makes it the missing
// middle term: deviation from the ODE is fluctuations, deviation from the 2D lattice
// run is spatial correlation.
//
// Sites are unlabelled without a lattice, so the state is just two pools of masses
// (positive and negative magnitudes). Picking among occupied sites rather than all L
// only rescales time.
//
// Usage: ./wellMixed <N> <sweeps> <q> <seed> [outDir]

#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <sstream>
#include <filesystem>

constexpr long long DEFAULT_N = 200000;
constexpr long long DEFAULT_SWEEPS = 1000;
constexpr double DEFAULT_Q = 0.8;
constexpr int RECORD_INTERVAL = 10;

std::mt19937_64 gen;

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

inline void removeAt(std::vector<long long> &v, size_t i)
{
    v[i] = v.back();
    v.pop_back();
}

inline size_t pick(size_t n)
{
    std::uniform_int_distribution<size_t> d(0, n - 1);
    return d(gen);
}

// One collision: a spot uniformly at random, then a partner of the same sign with
// probability q and of the opposite sign otherwise. The pair is drawn ~ n_i n_j
// either way, which is the constant kernel of the rate equation.
void update(std::vector<long long> &pos, std::vector<long long> &neg,
            long long N, double q, bool record, Hist &emis)
{
    size_t np = pos.size(), nn = neg.size();
    if (np == 0 || nn == 0) return;
    std::uniform_real_distribution<double> dr(0.0, 1.0);

    bool aPos = pick(np + nn) < np;
    std::vector<long long> &A = aPos ? pos : neg;
    size_t ia = pick(A.size());

    if (dr(gen) < q)                        // same sign: coalesce, i,j -> i+j
    {
        if (A.size() < 2) return;
        size_t ib;
        do { ib = pick(A.size()); } while (ib == ia);
        A[ia] += A[ib];                     // merge first, then swap-pop the partner
        removeAt(A, ib);
    }
    else                                    // opposite sign: annihilate
    {
        std::vector<long long> &B = aPos ? neg : pos;
        size_t ib = pick(B.size());
        long long a = A[ia], b = B[ib];
        if (record) emis.add(std::min(a, b));
        if (a > b)      { A[ia] = a - b; removeAt(B, ib); }
        else if (a < b) { B[ib] = b - a; removeAt(A, ia); }
        else            { removeAt(A, ia); removeAt(B, ib); }
    }

    // hold the spot count fixed with charge-neutral monomer pairs
    while ((long long)(pos.size() + neg.size()) < N)
    {
        pos.push_back(1);
        neg.push_back(1);
    }
}

int main(int argc, char *argv[])
{
    long long N = DEFAULT_N, sweeps = DEFAULT_SWEEPS;
    double q = DEFAULT_Q;
    unsigned seed = 1;
    std::string outDir = "outputs";
    if (argc > 1) N = std::stoll(argv[1]);
    if (argc > 2) sweeps = std::stoll(argv[2]);
    if (argc > 3) q = std::stod(argv[3]);
    if (argc > 4) seed = (unsigned)std::stoul(argv[4]);
    if (argc > 5) outDir = argv[5];
    gen.seed(seed);

    std::vector<long long> pos(N / 2, 1), neg(N / 2, 1);
    Hist spot, emis;
    long long rec = sweeps / 2;             // discard the first half as transient
    for (long long s = 0; s < sweeps; ++s)
    {
        bool record = s >= rec;
        for (long long i = 0; i < N; ++i) update(pos, neg, N, q, record, emis);
        if (record && s % RECORD_INTERVAL == 0)
        {
            for (long long v : pos) spot.add(v);
            for (long long v : neg) spot.add(v);
        }
    }

    std::filesystem::create_directories(outDir);
    std::ostringstream tag;
    tag << "N_" << N << "_q_" << q << "_seed_" << seed;
    std::ofstream sf(outDir + "/wmSpotSize_" + tag.str() + ".tsv");
    std::ofstream ef(outDir + "/wmEmission_" + tag.str() + ".tsv");
    spot.write(sf);
    emis.write(ef);
    return 0;
}
