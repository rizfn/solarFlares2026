import matplotlib.pyplot as plt
import numpy as np

def power_law(x, tau, alpha):
    return alpha * x**-tau

def plot_Neighbour():
    L = 64
    density = 0.2
    steps = 1000000

    tau_spot = 1.5
    tau_emission = 1.9

    spot_size_file_path = f'src/PUKPaper/outputs/2DNeighbour/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'
    emission_file_path = spot_size_file_path.replace('spotSize', 'emission')

    # Initialize an empty list to store spot sizes
    spot_sizes = []

    # Read the spot size file line by line
    with open(spot_size_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                sizes = list(map(int, filter(None, parts[1].split(','))))
                spot_sizes.extend(sizes)

    # Read the emission histogram file: two columns, "size\tcount"
    emission_data = np.loadtxt(emission_file_path, dtype=np.int64, ndmin=2)
    emission_sizes = emission_data[:, 0]
    emission_counts = emission_data[:, 1]

    # Convert the list to a numpy array
    spot_sizes = np.abs(np.array(spot_sizes))

    # Filter out non-positive values  todo: is this needed??
    spot_sizes = spot_sizes[spot_sizes > 0]
    emission_mask = emission_sizes > 0
    emission_sizes = emission_sizes[emission_mask]
    emission_counts = emission_counts[emission_mask]

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
    hist_emission, bin_edges_emission = np.histogram(emission_sizes, bins=bins_emission, weights=emission_counts)
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
    ax1.grid()
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
    ax2.grid()
    ax2.set_xlabel('Emission Size')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Emission Sizes')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(f'src/PUKPaper/plots/2DNeighbour/L_{L}_density_{density}_steps_{steps}.png', dpi=300)
    plt.show()


def plot_random():
    L = 128
    density = 0.2
    steps = 1000000

    tau_spot = 2
    tau_emission = 3

    spot_size_file_path = f'src/PUKPaper/outputs/2DRandom/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'  # Replace with your actual file path
    emission_file_path = spot_size_file_path.replace('spotSize', 'emission')

    # Initialize an empty list to store spot sizes
    spot_sizes = []

    # Read the spot size file line by line
    with open(spot_size_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                sizes = list(map(int, filter(None, parts[1].split(','))))
                spot_sizes.extend(sizes)

    # Read the emission histogram file: two columns, "size\tcount"
    emission_data = np.loadtxt(emission_file_path, dtype=np.int64, ndmin=2)
    emission_sizes = emission_data[:, 0]
    emission_counts = emission_data[:, 1]

    # Convert the list to a numpy array
    spot_sizes = np.abs(np.array(spot_sizes))

    # Filter out non-positive values
    spot_sizes = spot_sizes[spot_sizes > 0]
    emission_mask = emission_sizes > 0
    emission_sizes = emission_sizes[emission_mask]
    emission_counts = emission_counts[emission_mask]

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
    hist_emission, bin_edges_emission = np.histogram(emission_sizes, bins=bins_emission, weights=emission_counts)
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
    ax1.grid()
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
    ax2.grid()
    ax2.set_xlabel('Emission Size')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Emission Sizes')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(f'src/PUKPaper/plots/2DRandom/L_{L}_density_{density}_steps_{steps}.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    plot_Neighbour()
    # plot_random()