import matplotlib.pyplot as plt
import numpy as np

def load_spot_sizes(file_path):
    """Load spot sizes and timesteps from the given file."""
    timesteps = []
    spot_sizes = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                timesteps.append(int(parts[0]))  # First value is the timestep
                sizes = list(map(int, filter(None, parts[1].split(','))))
                spot_sizes.append(sizes)
    return timesteps, spot_sizes

def plot_max_min_over_time():
    # Parameters
    L = 128
    density = 0.2
    steps = 1000000

    # File paths
    neighbour_file_path = f'PUK_paper/outputs/2DNeighbour/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'
    random_file_path = f'PUK_paper/outputs/2DRandom/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'

    # 1D!!
    # L = 1024
    # neighbour_file_path = f'PUK_paper/outputs/1DNeighbour/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'
    # random_file_path = f'PUK_paper/outputs/1DRandom/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'


    # Load spot sizes and timesteps for both cases
    neighbour_timesteps, neighbour_spot_sizes = load_spot_sizes(neighbour_file_path)
    random_timesteps, random_spot_sizes = load_spot_sizes(random_file_path)

    # Compute max and min spot sizes over time
    neighbour_max_sizes = [max(sizes) if sizes else 0 for sizes in neighbour_spot_sizes]
    neighbour_min_sizes = [min(sizes) if sizes else 0 for sizes in neighbour_spot_sizes]

    random_max_sizes = [max(sizes) if sizes else 0 for sizes in random_spot_sizes]
    random_min_sizes = [min(sizes) if sizes else 0 for sizes in random_spot_sizes]

    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Plot max and min for the neighbour model
    ax1.plot(neighbour_timesteps, neighbour_max_sizes, label='Neighbour Max', color='blue')
    ax1.plot(neighbour_timesteps, neighbour_min_sizes, label='Neighbour Min', color='red')
    ax1.set_xlabel('Timesteps')
    ax1.set_ylabel('Spot Size')
    ax1.set_title('Neighbour Model: Max and Min Spot Sizes Over Time')
    ax1.legend()
    ax1.grid()

    # Plot max and min for the random model
    ax2.plot(random_timesteps, random_max_sizes, label='Random Max', color='green')
    ax2.plot(random_timesteps, random_min_sizes, label='Random Min', color='orange')
    ax2.set_xlabel('Timesteps')
    ax2.set_ylabel('Spot Size')
    ax2.set_title('Random Model: Max and Min Spot Sizes Over Time')
    ax2.legend()
    ax2.grid()

    # Adjust layout and save the plot
    plt.tight_layout()
    # plt.savefig(f'PUK_paper/plots/2DMaxMin/L_{L}_density_{density}_steps_{steps}.png', dpi=300)
    plt.show()



def plot_max_min_and_return_times():
    """Plot deviations and return time histograms for neighbour and random models."""
    # Parameters
    L = 128
    density = 0.2
    steps = 1000000

    # File paths
    neighbour_file_path = f'PUK_paper/outputs/2DNeighbour/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'
    random_file_path = f'PUK_paper/outputs/2DRandom/spotSize_L_{L}_density_{density}_steps_{steps}.tsv'

    # Load spot sizes and timesteps for both cases
    neighbour_timesteps, neighbour_spot_sizes = load_spot_sizes(neighbour_file_path)
    random_timesteps, random_spot_sizes = load_spot_sizes(random_file_path)

    # Helper function to calculate deviations and return times
    def calculate_deviation_and_return_times(timesteps, spot_sizes):
        max_abs_sizes = [max(abs(size) for size in sizes) if sizes else 0 for sizes in spot_sizes]
        avg_max_abs_size = np.mean(max_abs_sizes)
        deviations = [max_abs - avg_max_abs_size for max_abs in max_abs_sizes]
        # Calculate return times
        return_times = []
        current_return_time = 0
        for i in range(1, len(deviations)):
            current_return_time += 1
            if deviations[i] * deviations[i - 1] < 0:
                return_times.append(current_return_time)
                current_return_time = 0
        return deviations, return_times

    # Calculate deviations and return times for both cases
    neighbour_deviations, neighbour_return_times = calculate_deviation_and_return_times(neighbour_timesteps, neighbour_spot_sizes)
    random_deviations, random_return_times = calculate_deviation_and_return_times(random_timesteps, random_spot_sizes)

    # Create subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Top row: Deviations
    ax1, ax2 = axes[0]
    ax1.plot(neighbour_timesteps, neighbour_deviations, label='Deviation', color='blue')
    ax1.axhline(0, linestyle='--', color='black', label='Zero Line')
    ax1.set_xlabel('Timesteps')
    ax1.set_ylabel('Deviation')
    ax1.set_title('Neighbour Model: Deviation from Mean Max Abs Spot Size')
    ax1.legend()
    ax1.grid()

    ax2.plot(random_timesteps, random_deviations, label='Deviation', color='green')
    ax2.axhline(0, linestyle='--', color='black', label='Zero Line')
    ax2.set_xlabel('Timesteps')
    ax2.set_ylabel('Deviation')
    ax2.set_title('Random Model: Deviation from Mean Max Abs Spot Size')
    ax2.legend()
    ax2.grid()

    # Bottom row: Histograms
    ax3, ax4 = axes[1]
    bins = np.histogram_bin_edges(neighbour_return_times + random_return_times, bins=20)  # Shared bins    

    ax3.hist(neighbour_return_times, bins=bins, color='gray', edgecolor='black', alpha=0.7)
    ax3.set_xlabel('Return Time (Timesteps)')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Neighbour Model: Return Time Distribution')
    ax3.grid()

    ax4.hist(random_return_times, bins=bins, color='gray', edgecolor='black', alpha=0.7)
    ax4.set_xlabel('Return Time (Timesteps)')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Random Model: Return Time Distribution')
    ax4.grid()

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_max_min_over_time()
    plot_max_min_and_return_times()