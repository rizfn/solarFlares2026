// Instrumented copy of probabilisticNeighbour.cpp that follows the identity of every
// spot, so the coagulation cascade can be reconstructed after the fact.
//
// Each spot gets a unique id when it is injected, together with the site it was born on.
// When two spots meet, the survivor keeps the id of the heavier one and the other id is
// retired; the event is logged. From that log you can build
//   - the merger tree of any final spot (who ate whom, and when), and
//   - its drainage basin: the birth sites of every ancestor that fed into it.
//
// Written for short runs started from a random state: the whole point is to watch the
// cascade build up, and the event log is one line per collision, so it is not meant to
// run for 1e5 sweeps.
//
// usage: ./lineage L rho steps p seed outDir
#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <sstream>
#include <filesystem>

std::random_device rd;
std::mt19937 gen(rd());

long long nextId = 0;

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

// give a freshly injected spot an identity and record where it entered the system
void birth(int loc, int L, std::vector<long long> &id, std::ofstream &ev, long long step)
{
    id[loc] = nextId++;
    ev << "B\t" << step << "\t" << id[loc] << "\t" << (loc % L) << "\t" << (loc / L) << "\n";
}

void addNeighbourPair(std::vector<long long> &s, int L, std::vector<int> &f,
                      std::vector<int> &pos, std::vector<long long> &id,
                      std::ofstream &ev, long long step)
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
                if (s[nb] == 0) { s[nb] = 1; addFilled(nb, f, pos); birth(nb, L, id, ev, step); posF = true; break; }
        }
        else if (s[loc] < 0 && !negF)
        {
            for (int nb : neighbours(loc, L))
                if (s[nb] == 0) { s[nb] = -1; addFilled(nb, f, pos); birth(nb, L, id, ev, step); negF = true; break; }
        }
    }
    std::uniform_int_distribution<> dl(0, L * L - 1);
    if (!posF) { int r; do { r = dl(gen); } while (s[r] != 0); s[r] = 1; addFilled(r, f, pos); birth(r, L, id, ev, step); }
    if (!negF) { int r; do { r = dl(gen); } while (s[r] != 0); s[r] = -1; addFilled(r, f, pos); birth(r, L, id, ev, step); }
}

void addRandomPair(std::vector<long long> &s, int L, std::vector<int> &f,
                   std::vector<int> &pos, std::vector<long long> &id,
                   std::ofstream &ev, long long step)
{
    std::uniform_int_distribution<> dl(0, L * L - 1);
    int a; do { a = dl(gen); } while (s[a] != 0);
    s[a] = 1; addFilled(a, f, pos); birth(a, L, id, ev, step);
    int b; do { b = dl(gen); } while (s[b] != 0);
    s[b] = -1; addFilled(b, f, pos); birth(b, L, id, ev, step);
}

void update(std::vector<long long> &s, int L, int N, double p, std::vector<int> &f,
            std::vector<int> &pos, std::vector<long long> &id,
            std::ofstream &ev, long long step)
{
    std::uniform_int_distribution<> di(0, (int)f.size() - 1);
    int idx = di(gen), loc = f[idx];
    long long val = s[loc], vid = id[loc];
    std::uniform_real_distribution<> dr(0.0, 1.0);
    int x = loc % L, y = loc / L, dst;
    double r = dr(gen);
    if (r < 0.25) dst = ((y - 1 + L) % L) * L + x;
    else if (r < 0.5) dst = ((y + 1) % L) * L + x;
    else if (r < 0.75) dst = y * L + (x - 1 + L) % L;
    else dst = y * L + (x + 1) % L;

    long long dv = s[dst], did = id[dst];
    removeAt(idx, f, pos);
    if (dv != 0)
    {
        long long res = val + dv;
        if (val * dv > 0)
            // coagulation: the target keeps its identity, the mover is absorbed
            ev << "C\t" << step << "\t" << did << "\t" << vid << "\t"
               << std::llabs(dv) << "\t" << std::llabs(val) << "\n";
        else
            ev << "A\t" << step << "\t" << did << "\t" << vid << "\t"
               << std::llabs(dv) << "\t" << std::llabs(val) << "\n";
        if (res == 0) { removeLoc(dst, f, pos); id[dst] = -1; }
        // the survivor is whichever spot was heavier; a lighter target loses its identity
        else if (std::llabs(val) > std::llabs(dv)) id[dst] = vid;
    }
    else { addFilled(dst, f, pos); id[dst] = vid; }
    s[dst] += val; s[loc] -= val;
    if (s[loc] == 0) id[loc] = -1;

    if ((int)f.size() < N)
    {
        if (dr(gen) < p) addNeighbourPair(s, L, f, pos, id, ev, step);
        else addRandomPair(s, L, f, pos, id, ev, step);
    }
}

int main(int argc, char *argv[])
{
    long long L = 64, steps = 1000;
    double rho = 0.2, p = 1.0;
    unsigned seed = 1;
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
    std::ofstream ev(outDir + "/lineage_" + tag.str() + ".tsv");
    std::ofstream fin(outDir + "/lineageFinal_" + tag.str() + ".tsv");
    ev << "# B step id x y | C/A step survivorId moverId mTarget mMover\n";

    int N = (int)((L * L * rho) / 2) * 2;
    std::vector<long long> s(L * L, 0), id(L * L, -1);
    std::fill(s.begin(), s.begin() + N / 2, 1);
    std::fill(s.begin() + N / 2, s.begin() + N, -1);
    std::shuffle(s.begin(), s.end(), gen);

    std::vector<int> f; f.reserve(2 * N + 8);
    std::vector<int> pos(L * L, -1);
    for (int i = 0; i < L * L; ++i)
        if (s[i]) { addFilled(i, f, pos); id[i] = nextId++;
                    ev << "B\t0\t" << id[i] << "\t" << (i % L) << "\t" << (i / L) << "\n"; }

    for (long long step = 0; step < steps; ++step)
        for (int i = 0; i < L * L; ++i) update(s, L, N, p, f, pos, id, ev, step);

    // final state: every surviving spot with its identity, so the tree has a root
    fin << "# id\tx\ty\tmass\n";
    for (int loc : f)
        fin << id[loc] << "\t" << (loc % L) << "\t" << (loc / L) << "\t" << s[loc] << "\n";
    return 0;
}
