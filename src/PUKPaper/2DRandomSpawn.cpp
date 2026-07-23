#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <iomanip>
#include <filesystem>

constexpr long long DEFAULT_L = 64;
constexpr double DEFAULT_DENSITY = 0.2;
constexpr long long DEFAULT_STEPS_PER_LATTICEPOINT = 1000000;
constexpr int RECORDING_INTERVAL = 100;

std::random_device rd;
std::mt19937 gen(rd());

struct EmissionHist
{
    static constexpr long long LIMIT = 1LL << 22;
    std::vector<long long> small;
    std::unordered_map<long long, long long> large;

    EmissionHist() : small(1024, 0) {}

    void add(long long e)
    {
        if (e >= LIMIT)
        {
            ++large[e];
            return;
        }
        if (e >= static_cast<long long>(small.size()))
        {
            long long newSize = static_cast<long long>(small.size());
            while (newSize <= e)
                newSize *= 2;
            small.resize(std::min(newSize, LIMIT), 0);
        }
        ++small[e];
    }

    void write(std::ofstream &file) const
    {
        file << "# emissionSize\tcount\n";
        for (long long i = 0; i < static_cast<long long>(small.size()); ++i)
            if (small[i] != 0)
                file << i << "\t" << small[i] << "\n";
        std::vector<long long> sizes;
        for (const auto &entry : large)
            sizes.push_back(entry.first);
        std::sort(sizes.begin(), sizes.end());
        for (long long size : sizes)
            file << size << "\t" << large.at(size) << "\n";
    }
};

void addFilled(int loc, std::vector<int> &filledLocs, std::vector<int> &pos)
{
    pos[loc] = static_cast<int>(filledLocs.size());
    filledLocs.push_back(loc);
}

void removeFilledAt(int idx, std::vector<int> &filledLocs, std::vector<int> &pos)
{
    int loc = filledLocs[idx];
    int last = filledLocs.back();
    filledLocs[idx] = last;
    pos[last] = idx;
    filledLocs.pop_back();
    pos[loc] = -1;
}

void removeFilledLoc(int loc, std::vector<int> &filledLocs, std::vector<int> &pos)
{
    removeFilledAt(pos[loc], filledLocs, pos);
}

void addNewSpots(std::vector<long long> &state, int L, int N, std::vector<int> &filledLocs, std::vector<int> &pos)
{
    bool posFound = false, negFound = false;
    std::uniform_int_distribution<> dis(0, L * L - 1);

    while (!posFound || !negFound)
    {
        int randomLoc = dis(gen);
        if (state[randomLoc] == 0)
        {
            if (!posFound)
            {
                state[randomLoc] = 1;
                addFilled(randomLoc, filledLocs, pos);
                posFound = true;
            }
            else if (!negFound)
            {
                state[randomLoc] = -1;
                addFilled(randomLoc, filledLocs, pos);
                negFound = true;
            }
        }
    }
}

int update(std::vector<long long> &state, int L, int N, EmissionHist &emissionHist, std::vector<int> &filledLocs, std::vector<int> &pos)
{
    std::uniform_int_distribution<> dis(0, filledLocs.size() - 1);
    int idx = dis(gen);
    int spotLoc = filledLocs[idx];
    long long stateVal = state[spotLoc];
    std::uniform_real_distribution<> disReal(0.0, 1.0);
    int newSpotLoc;
    double r = disReal(gen);
    if (r < 0.25)
        newSpotLoc = ((spotLoc / L - 1 + L) % L) * L + spotLoc % L; // up
    else if (r < 0.5)
        newSpotLoc = ((spotLoc / L + 1) % L) * L + spotLoc % L;     // down
    else if (r < 0.75)
        newSpotLoc = spotLoc / L * L + (spotLoc % L - 1 + L) % L;   // left
    else
        newSpotLoc = spotLoc / L * L + (spotLoc % L + 1) % L;       // right

    long long currentNewLocVal = state[newSpotLoc];
    long long emiss = 0;
    removeFilledAt(idx, filledLocs, pos);
    if (currentNewLocVal != 0)
    {
        if (stateVal * currentNewLocVal < 0)
        {
            emiss = std::min(std::abs(stateVal), std::abs(currentNewLocVal));
            emissionHist.add(emiss);
        }
        if (stateVal == -currentNewLocVal)
        {
            removeFilledLoc(newSpotLoc, filledLocs, pos);
        }
    }
    else
    {
        addFilled(newSpotLoc, filledLocs, pos);
    }
    state[newSpotLoc] += stateVal;
    state[spotLoc] -= stateVal;

    if (static_cast<int>(filledLocs.size()) < N)
    {
        addNewSpots(state, L, N, filledLocs, pos);
    }

    return emiss;
}

void run(std::ofstream &spotSizeFile, std::ofstream &emissionFile, int L, int N, int stepsPerLatticepoint, int recordingStep)
{
    std::vector<long long> state(L * L, 0);
    std::fill(state.begin(), state.begin() + N / 2, 1);
    std::fill(state.begin() + N / 2, state.begin() + N, -1);
    std::shuffle(state.begin(), state.end(), gen);

    std::vector<int> filledLocs;
    filledLocs.reserve(2 * N + 8);
    std::vector<int> pos(L * L, -1);
    for (int i = 0; i < L * L; ++i)
    {
        if (state[i] != 0)
            addFilled(i, filledLocs, pos);
    }

    EmissionHist emissionHist;

    for (int step = 0; step < stepsPerLatticepoint; ++step)
    {
        for (int i = 0; i < L * L; ++i)
        {
            update(state, L, N, emissionHist, filledLocs, pos);
        }
        if ((step >= recordingStep) && (step % RECORDING_INTERVAL == 0))
        {
            spotSizeFile << step << "\t";
            for (long long val : state)
            {
                if (val != 0)
                {
                    spotSizeFile << val << ",";
                }
            }
            spotSizeFile.seekp(-1, std::ios_base::cur); // Remove the last comma
            spotSizeFile << "\n";
        }
        if (step % 1000 == 0)
            std::cout << "Progress: " << std::fixed << std::setprecision(2) << static_cast<double>(step) / stepsPerLatticepoint * 100 << "%\r" << std::flush;
    }

    emissionHist.write(emissionFile);
    emissionFile.close();
    spotSizeFile.close();
}

int main(int argc, char *argv[])
{
    long long L = DEFAULT_L;
    double density = DEFAULT_DENSITY;
    long long stepsPerLatticepoint = DEFAULT_STEPS_PER_LATTICEPOINT;

    if (argc > 1)
        L = std::stoll(argv[1]);
    if (argc > 2)
        density = std::stod(argv[2]);
    if (argc > 3)
        stepsPerLatticepoint = std::stoll(argv[3]);

    long long recordingStep = stepsPerLatticepoint / 2; // start recording after half the steps
    int N = static_cast<int>((L * L * density) / 2) * 2;

    std::string exePath = argv[0];
    std::string exeDir = std::filesystem::path(exePath).parent_path().string();
    std::filesystem::create_directories(exeDir + "/outputs/2DRandom");
    std::ostringstream spotSizePathStream;
    spotSizePathStream << exeDir << "/outputs/2DRandom/spotSize_L_" << L << "_density_" << density << "_steps_" << stepsPerLatticepoint << ".tsv";
    std::string spotSizePath = spotSizePathStream.str();
    std::ofstream spotSizeFile;
    spotSizeFile.open(spotSizePath);

    std::string emissionPath = spotSizePath.replace(spotSizePath.find("spotSize"), 8, "emission");
    std::ofstream emissionFile;
    emissionFile.open(emissionPath);

    run(spotSizeFile, emissionFile, L, N, stepsPerLatticepoint, recordingStep);

    return 0;
}
