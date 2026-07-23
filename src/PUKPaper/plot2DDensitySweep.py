import matplotlib.pyplot as plt
import numpy as np

def power_law(x, tau, alpha):
    return alpha * x**-tau

def read_emission_histogram(file_path):
    # Read the emission histogram file: two columns, "size\tcount"
    data = np.loadtxt(file_path, dtype=np.int64, ndmin=2)
    sizes = data[:, 0]
    counts = data[:, 1]
    mask = sizes > 0
    return sizes[mask], counts[mask]

def fit_tau(sizes, counts, xmin):
    # Maximum likelihood exponent for P(s) ~ s^-tau, s >= xmin
    mask = sizes >= xmin
    sizes = sizes[mask]
    counts = counts[mask]
    return 1 + counts.sum() / np.sum(counts * np.log(sizes / xmin))

def plot_Neighbour():
    L = 128
    steps = 1000000
    densities = [0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9]
    xmin = 10

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(densities)))

    for density, color in zip(densities, colors):
        emission_file_path = f'src/PUKPaper/outputs/2DNeighbour/emission_L_{L}_density_{density}_steps_{steps}.tsv'
        emission_sizes, emission_counts = read_emission_histogram(emission_file_path)
        tau = fit_tau(emission_sizes, emission_counts, xmin)

        # Create log-spaced bins for emission sizes
        bins_emission = np.geomspace(np.min(emission_sizes), np.max(emission_sizes), num=50)

        # Compute the histogram for emission sizes
        hist_emission, bin_edges_emission = np.histogram(emission_sizes, bins=bins_emission, weights=emission_counts)
        bin_widths_emission = np.diff(bin_edges_emission)
        hist_normalized_emission = hist_emission / bin_widths_emission

        ax.plot(bin_edges_emission[:-1], hist_normalized_emission, marker='x', linestyle='none',
                color=color, label=f'$\\rho$={density}, $\\tau$={tau:.2f}')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid()
    ax.set_xlabel('Emission Size')
    ax.set_ylabel('Frequency')
    ax.set_title('Neighbour Model: Emission Sizes vs Density')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'src/PUKPaper/plots/2DNeighbour/emissionPDF_densitySweep_L_{L}_steps_{steps}.png', dpi=300)
    plt.show()


def plot_random():
    L = 128
    steps = 1000000
    densities = [0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9]
    xmin = 10

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(densities)))

    for density, color in zip(densities, colors):
        emission_file_path = f'src/PUKPaper/outputs/2DRandom/emission_L_{L}_density_{density}_steps_{steps}.tsv'
        emission_sizes, emission_counts = read_emission_histogram(emission_file_path)
        tau = fit_tau(emission_sizes, emission_counts, xmin)

        # Create log-spaced bins for emission sizes
        bins_emission = np.geomspace(np.min(emission_sizes), np.max(emission_sizes), num=50)

        # Compute the histogram for emission sizes
        hist_emission, bin_edges_emission = np.histogram(emission_sizes, bins=bins_emission, weights=emission_counts)
        bin_widths_emission = np.diff(bin_edges_emission)
        hist_normalized_emission = hist_emission / bin_widths_emission

        ax.plot(bin_edges_emission[:-1], hist_normalized_emission, marker='x', linestyle='none',
                color=color, label=f'$\\rho$={density}, $\\tau$={tau:.2f}')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid()
    ax.set_xlabel('Emission Size')
    ax.set_ylabel('Frequency')
    ax.set_title('Random Model: Emission Sizes vs Density')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'src/PUKPaper/plots/2DRandom/emissionPDF_densitySweep_L_{L}_steps_{steps}.png', dpi=300)
    plt.show()


def plot_tau_vs_density():
    L = 128
    steps = 1000000
    densities = [0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9]
    xmin = 10

    fig, ax = plt.subplots(figsize=(8, 6))

    for model, label, marker in [('2DNeighbour', 'Neighbour Model', 'o'), ('2DRandom', 'Random Model', 's')]:
        taus = []
        for density in densities:
            emission_file_path = f'src/PUKPaper/outputs/{model}/emission_L_{L}_density_{density}_steps_{steps}.tsv'
            emission_sizes, emission_counts = read_emission_histogram(emission_file_path)
            taus.append(fit_tau(emission_sizes, emission_counts, xmin))
        ax.plot(densities, taus, marker=marker, label=label)

    ax.grid()
    ax.set_xlabel('Density')
    ax.set_ylabel('Emission Exponent $\\tau$')
    ax.set_title('Emission Exponent vs Density')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'src/PUKPaper/plots/2DDensitySweep/tau_vs_density_L_{L}_steps_{steps}.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    plot_Neighbour()
    # plot_random()
    # plot_tau_vs_density()
