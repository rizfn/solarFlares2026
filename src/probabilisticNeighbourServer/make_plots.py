# Figures from the server sweep. Reuses the analysis in ../probabilisticNeighbour/plots.py
# with the input and output directories pointed here.
import os
import sys

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(DIR, "..", "probabilisticNeighbour"))

import plots as P

P.OUT = os.path.join(DIR, "outputs")
P.PLOTS = os.path.join(DIR, "plots")

LS = (128, 256, 512, 1024)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"

    if which in ("all", "cheap"):
        os.makedirs(os.path.join(P.PLOTS, "snapshots"), exist_ok=True)
        os.makedirs(os.path.join(P.PLOTS, "histograms"), exist_ok=True)
        os.makedirs(os.path.join(P.PLOTS, "exponents"), exist_ok=True)
        os.makedirs(os.path.join(P.PLOTS, "finiteSize"), exist_ok=True)
        P.plot_histograms(L=128, rhos=(0.2, 0.4, 0.6, 0.8), ps=(0.0, 0.5, 1.0))
        P.plot_histograms(L=512, rhos=(0.2,), ps=(0.0, 0.5, 1.0))
        P.plot_finite_size(rho=0.2, p=1.0, Ls=LS)
        P.plot_snapshots(L=512, rhos=(0.2,), ps=(0.4, 0.6, 0.8))
        P.plot_snapshots(L=128, rhos=(0.2, 0.6), ps=(0.0, 0.5, 1.0))
        P.plot_max_trajectories(L=512, rho=0.2)
        P.plot_max_finite_size(rho=0.2, Ls=LS, ps=(0.0, 0.5, 1.0))
        P.plot_correlation_function(rho=0.2, Ls=(256, 512, 1024), ps=(0.6, 0.8, 0.9, 1.0))

    if which in ("all", "voronoi"):
        # the expensive ones: a Delaunay per snapshot, ~5 s at L=512 and ~20 s at L=1024
        P.plot_exponents_vs_p(L=128, rho=0.2)
        P.plot_exponents_vs_p(L=512, rho=0.2)
        P.plot_voronoi_same_sign(L=128, rhos=(0.2, 0.4, 0.6, 0.8))
        P.plot_snapshots_voronoi(L=128, rhos=(0.2, 0.6), ps=(0.0, 0.5, 1.0))
        P.plot_qv_finite_size(rho=0.2, Ls=LS)
    print("done:", which)
