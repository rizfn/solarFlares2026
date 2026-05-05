import matplotlib.pyplot as plt
import numpy as np
import sys
import os

plt.style.use('seaborn-whitegrid')
plt.rcParams['figure.dpi'] = 300

def power_law(x, tau, alpha):
    return alpha * x**-tau

def main():
    L = 1024
    density = 0.2
    steps = 1000000

    tau_spot = 1.3
    tau_emission = 1.6
    tau_waiting_time = 1

    spot_size_file_path = f'PUK_paper/outputs/1DNeighbour/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'  # Replace with your actual file path
    # spot_size_file_path = f'PUK_paper/outputs/1DRandom/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'  # Replace with your actual file path
    emission_file_path = spot_size_file_path.replace('spotSize', 'emission')

    # Initialize empty lists to store spot sizes and emission sizes
    spot_sizes = []
    emission_sizes = []
    emission_steps = []

    # Read the spot size file line by line
    with open(spot_size_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                sizes = list(map(int, filter(None, parts[1].split(','))))
                spot_sizes.extend(sizes)

    # Read the emission file line by line
    with open(emission_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                emission_steps.append(float(parts[0]))
                emission_sizes.append(int(parts[1]))

    # Convert the lists to numpy arrays
    spot_sizes = np.abs(np.array(spot_sizes))
    emission_sizes = np.abs(np.array(emission_sizes))
    emission_steps = np.array(emission_steps)

    # Filter out non-positive values
    spot_sizes = spot_sizes[spot_sizes > 0]
    emission_sizes = emission_sizes[emission_sizes > 0]

    # Create log-spaced bins for spot sizes
    min_size_spot = np.min(spot_sizes)
    max_size_spot = np.max(spot_sizes)
    bins_spot = np.geomspace(min_size_spot, max_size_spot, num=50)

    # Compute the histogram for spot sizes
    hist_spot, bin_edges_spot = np.histogram(spot_sizes, bins=bins_spot)
    bin_widths_spot = np.diff(bin_edges_spot)
    hist_normalized_spot = hist_spot / bin_widths_spot

    # Create log-spaced bins for emission sizes
    min_size_emission = np.min(emission_sizes)
    max_size_emission = np.max(emission_sizes)
    bins_emission = np.geomspace(min_size_emission, max_size_emission, num=50)

    # Compute the histogram for emission sizes
    hist_emission, bin_edges_emission = np.histogram(emission_sizes, bins=bins_emission)
    bin_widths_emission = np.diff(bin_edges_emission)
    hist_normalized_emission = hist_emission / bin_widths_emission

    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Plot the histogram for spot sizes
    ax1.plot(bin_edges_spot[:-1], hist_normalized_spot, marker='x', linestyle='none', label='Histogram')
    x_spot = np.linspace(min_size_spot, max_size_spot, 100)
    y_spot = power_law(x_spot, tau_spot, alpha=hist_normalized_spot.max())
    ax1.plot(x_spot, y_spot, label=f'$\\tau$={tau_spot} Power Law', linestyle='--')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Spot Size')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Spot Sizes')
    ax1.legend()

    # Plot the histogram for emission sizes
    ax2.plot(bin_edges_emission[:-1], hist_normalized_emission, marker='x', linestyle='none', label='Histogram')
    x_emission = np.linspace(min_size_emission, max_size_emission, 100)
    y_emission = power_law(x_emission, tau_emission, alpha=hist_normalized_emission.max())
    ax2.plot(x_emission, y_emission, label=f'$\\tau$={tau_emission} Power Law', linestyle='--')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Emission Size')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Emission Sizes')
    ax2.legend()

    # # Calculate waiting times for different large emission sizes
    # large_emission_sizes = [10, 100, 1000]
    # waiting_times = {size: [] for size in large_emission_sizes}

    # for size in large_emission_sizes:
    #     large_emission_indices = np.where(emission_sizes >= size)[0]
    #     if len(large_emission_indices) > 1:
    #         waiting_times[size] = np.diff(emission_steps[large_emission_indices])

    # # Plot the waiting time distributions
    # for size in large_emission_sizes:
    #     if len(waiting_times[size]) > 0:
    #         positive_waiting_times = waiting_times[size][waiting_times[size] > 0]
    #         if len(positive_waiting_times) > 0:
    #             min_waiting_time = np.min(positive_waiting_times)
    #             max_waiting_time = np.max(positive_waiting_times)
    #             bins_waiting = np.linspace(min_waiting_time, max_waiting_time, num=50)
    #             hist_waiting, bin_edges_waiting = np.histogram(positive_waiting_times, bins=bins_waiting)
    #             bin_widths_waiting = np.diff(bin_edges_waiting)
    #             hist_normalized_waiting = hist_waiting / bin_widths_waiting
    #             ax3.plot(bin_edges_waiting[:-1], hist_normalized_waiting, marker='x', linestyle='none', label=f'Size >= {size}')
    # # ax3.set_xscale('log')
    # ax3.set_yscale('log')
    # ax3.set_xlabel('Waiting Time')
    # ax3.set_ylabel('Frequency')
    # ax3.set_title('Waiting Time Distributions')
    # ax3.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()