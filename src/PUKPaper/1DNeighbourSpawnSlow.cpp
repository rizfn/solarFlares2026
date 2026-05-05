#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <iomanip>
#include <filesystem>

constexpr long long DEFAULT_L = 1024;
constexpr double DEFAULT_DENSITY = 0.2;
constexpr long long DEFAULT_STEPS_PER_LATTICEPOINT = 1000000;
constexpr int RECORDING_INTERVAL = 100;

std::random_device rd;
std::mt19937 gen(rd());

std::vector<int> neighbours(int index, int L) {
    return {(index - 1 + L) % L, (index + 1) % L};
}

void addNewSpots(std::vector<int>& state, int L, int N) { // todo: don't use randomloc, use random from filledLocs
    bool posNbrFound = false, negNbrFound = false;
    std::uniform_int_distribution<> dis(0, L - 1);
    while (!posNbrFound || !negNbrFound) {
        int randomLoc = dis(gen);
        if (state[randomLoc] != 0) continue;
        for (int nbr : neighbours(randomLoc, L)) {
            if (!posNbrFound && state[nbr] > 0) {
                state[randomLoc] += 1;
                posNbrFound = true;
                break;
            }
            if (!negNbrFound && state[nbr] < 0) {
                state[randomLoc] -= 1;
                negNbrFound = true;
                break;
            }
        }
    }
}

int update(std::vector<int>& state, int L, int N, std::ofstream& emissionFile, double currentStep) {
    std::vector<int> filledLocs;  // todo: try not to recalculate and just update filledLocs
    for (int i = 0; i < L; ++i) {
        if (state[i] != 0) filledLocs.push_back(i);
    }
    int N_filled = filledLocs.size();
    std::uniform_int_distribution<> dis(0, N_filled - 1);
    int spotLoc = filledLocs[dis(gen)];
    int stateVal = state[spotLoc];
    std::uniform_real_distribution<> disReal(0.0, 1.0);
    int newSpotLoc = (disReal(gen) < 0.5) ? (spotLoc - 1 + L) % L : (spotLoc + 1) % L;
    int currentNewLocVal = state[newSpotLoc];
    int emiss = 0;
    if (currentNewLocVal != 0) {
        if (stateVal * currentNewLocVal < 0) {
            emiss = std::min(std::abs(stateVal), std::abs(currentNewLocVal));
            emissionFile << currentStep << "\t" << emiss << "\n";
        }
        if (stateVal == -currentNewLocVal) {
            N_filled -= 2;
        } else {
            N_filled -= 1;
        }
    }
    state[newSpotLoc] += stateVal;
    state[spotLoc] -= stateVal;

    if (N_filled < N - 1) {
        addNewSpots(state, L, N);
    }

    return emiss;
}

void run(std::ofstream& spotSizeFile, std::ofstream& emissionFile, int L, int N, int stepsPerLatticepoint, int recordingStep) {
    std::vector<int> state(L, 0);
    std::fill(state.begin(), state.begin() + N / 2, 1);
    std::fill(state.begin() + N / 2, state.begin() + N, -1);
    std::shuffle(state.begin(), state.end(), gen);

    for (int step = 0; step < stepsPerLatticepoint; ++step) {
        for (int i = 0; i < L; ++i) {
            double currentStep = step + static_cast<double>(i) / L;
            update(state, L, N, emissionFile, currentStep);
        }
        if (step % RECORDING_INTERVAL == 0) {
            spotSizeFile << step << "\t";
            for (int val : state) {
                if (val != 0) {
                    spotSizeFile << val << ",";
                }
            }
            spotSizeFile.seekp(-1, std::ios_base::cur); // Remove the last comma
            spotSizeFile << "\n";
        }
        std::cout << "Progress: " << std::fixed << std::setprecision(2) << static_cast<double>(step) / stepsPerLatticepoint * 100 << "%\r" << std::flush;
    }

    emissionFile.close();
    spotSizeFile.close();
}

int main(int argc, char *argv[]) {
    long long L = DEFAULT_L;
    double density = DEFAULT_DENSITY;
    long long stepsPerLatticepoint = DEFAULT_STEPS_PER_LATTICEPOINT;

    if (argc > 1)
        L = std::stoll(argv[1]);
    if (argc > 2)
        density = std::stod(argv[2]);
    if (argc > 3)
        stepsPerLatticepoint = std::stoll(argv[3]);

    long long recordingStep = stepsPerLatticepoint / 2;  // start recording after half the steps
    int N = static_cast<int>((L * density) / 2) * 2;

    std::string exePath = argv[0];
    std::string exeDir = std::filesystem::path(exePath).parent_path().string();
    std::ostringstream spotSizePathStream;
    spotSizePathStream << exeDir << "/outputs/1DNeighbour/spotSize_L_" << L << "_density_" << density << "_steps_" << stepsPerLatticepoint << ".tsv";
    std::string spotSizePath = spotSizePathStream.str();
    std::ofstream spotSizeFile;
    spotSizeFile.open(spotSizePath);

    std::string emissionPath = spotSizePath.replace(spotSizePath.find("spotSize"), 8, "emission");
    std::ofstream emissionFile;
    emissionFile.open(emissionPath);

    run(spotSizeFile, emissionFile, L, N, stepsPerLatticepoint, recordingStep);

    return 0;
}