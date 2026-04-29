
def plot_energies(energies, filename='energy_analysis'):

    import numpy as np
    import matplotlib.pyplot as plt

    n_cycles = len(energies)
    cycles = range(1, n_cycles+1)
    energy_avrg = np.mean(energies)
    x_ticks = np.arange(0, n_cycles+1, 10)

    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.linewidth": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "legend.fontsize": 9
    })

    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # Plot 1: Energy vs Cycle
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.axhline(y=energy_avrg, color='gray', linestyle=':', linewidth=1, alpha=0.9, label=f'Average energy: {energy_avrg:.2f}')
    ax1.plot(cycles, energies, linewidth=1.5, color='#2563eb')
    ax1.scatter(cycles, energies, color='#2563eb', s=10.0)
    #ax1.set_xticks(x_ticks, x_ticks)
    ax1.set_xlabel('Cycle')
    ax1.set_ylabel(r"$E$ [eV]")
    ax1.set_title('Energy Landscape', pad=10)
    #ax1.legend()

    # Plot 2: Sorted energies (global minima)
    sorted_energies = np.sort(energies)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.axhline(y=energy_avrg, color='gray', linestyle=':', linewidth=1, alpha=0.9, label=f'Average energy: {energy_avrg:.2f}')
    ax2.plot(cycles, sorted_energies, linewidth=1.5, color='#dc2626')
    ax2.scatter(cycles, sorted_energies, color='#dc2626', s=5.0)
    #ax2.set_xticks(x_ticks, x_ticks)
    ax2.set_xlabel('Rank')
    ax2.set_ylabel(r"$E$ [eV]")
    ax2.set_title('Sorted Energy Landscape', pad=10)
    ax2.legend()

    plt.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{filename}.pdf', bbox_inches='tight')


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = 'combined_ens.txt'

    import numpy as np
    from combine_energies import extract_energies
    energies = extract_energies(filename)
    energies = np.array(energies)

    plot_energies(energies, 'energy_analysis')