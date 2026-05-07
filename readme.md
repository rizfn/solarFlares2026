## Solar flares

This project aims at simplifying [Hughes et al's model](https://doi.org/10.1103/PhysRevLett.90.131101) for self-organaized criticality in coronary mass emission dynamics.


### Project Structure:

- `src` contains the source code
- `docs` contains interactive visualizations and presentations, if any exist


### Getting started

This repository uses C++ code for the simulations, and python for visualizations. 

The python packages are managed by [uv](https://docs.astral.sh/uv/). First, be sure that you have it installed! Once you do, run
```
uv venv
```
to create a virtual environment, followed by
```
uv sync
```
to automatically install all the required libraries.