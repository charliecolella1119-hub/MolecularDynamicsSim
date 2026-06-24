# Molecular Dynamics Simulator

A real-time C++ molecular dynamics sandbox for exploring Lennard-Jones particles, phase behavior, and statistical-mechanics observables. The project uses SFML for windowing and 2D visualization, with the simulation structured around periodic boundaries, pairwise Lennard-Jones forces, and velocity Verlet integration.


## Features

- Real-time Lennard-Jones particle simulation
- Velocity Verlet time integration
- Periodic boundary conditions with minimum-image distance
- Gas, liquid, crystal, hot gas, dense fluid, and melting-crystal presets
- Temperature, kinetic energy, potential energy, total energy, and pressure estimates
- Maxwell-Boltzmann velocity histogram overlay
- Normalized radial distribution function, `g(r)`
- Pressure vs temperature phase-diagram sampling
- Optional particle trails for motion inspection
- CSV export of thermodynamic time-series data

## Controls

| Key | Action |
| --- | --- |
| `1` | Gas preset |
| `2` | Liquid preset |
| `3` | Crystal preset |
| `4` | Hot gas experiment |
| `5` | Dense fluid experiment |
| `6` | Melting crystal experiment |
| `Up` / `Down` | Heat or cool the system |
| `T` | Temperature history graph |
| `G` | Radial distribution graph |
| `V` | Velocity histogram with Maxwell-Boltzmann fit |
| `P` | Pressure-temperature phase diagram |
| `H` | Toggle particle trails |
| `C` | Export `thermodynamic_log.csv` |
| `Space` | Pause or resume |
| `R` | Reset current preset |
| `Esc` | Quit |

## Build

Requirements:

- CMake 3.16 or newer
- C++17 compiler
- SFML 3

```bash
cmake -S . -B build
cmake --build build
./build/MolecularDynamicsSim
```

The app expects `assets/Arial.ttf` to be available relative to the build directory, as included in this repository.

## Scientific Notes

Particles interact through a Lennard-Jones potential,

```text
U(r) = 4 epsilon [ (sigma / r)^12 - (sigma / r)^6 ]
```

which captures a short-range repulsive core and longer-range attraction. The simulation uses periodic boundaries and the minimum-image convention so particles leaving one side of the box re-enter from the opposite side.

The displayed pressure is a two-dimensional virial-style estimate. The velocity histogram is compared against a two-dimensional Maxwell-Boltzmann speed distribution. The radial distribution function, `g(r)`, is normalized against the expected particle count in each annular shell, so a well-mixed gas tends toward `g(r) = 1` at larger distances.

