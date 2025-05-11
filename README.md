Three-Band Hubbard Model for Cuprates — Tight-Binding & Green’s Function Analysis

This Python project implements a three-band Hubbard model for cuprates using a tight-binding approach for U=0 (Hubbard parameter).
It describes a 2D lattice composed of copper-oxygen clusters (N-mers) and includes orbital degrees of freedom
(Cu 3d and O 2p orbitals: d, px, py) with various hopping terms: nearest-neighbor copper-oxygen hopping $\( t_{pd} \)$,
oxygen-oxygen hopping $\( t_{pp} \)$, and next-nearest-neighbor hopping $\( t_{pp'} \)$. 
The model also supports spatial modulation of the on-site energies to simulate structural or electronic ordering (e.g. CDWs), including,
yet not limited to, local effects on on-site copper energies (cosine-modulation) due to apical oxygen displacement.
The code is easily adaptable to check for other local effects on the hamiltonian parameters.

Key Features:
- Implements a realistic three-band model for high-$\( T_c \)$ cuprate superconductors
- Includes tight-binding parameters: $\( t_{pd} \)$, $\( t_{pp} \)$, $\( t_{pp'} \)$
- Allows spatial modulation of $\( \epsilon_d \)$ (on-site copper energy) via cosine profile
- Constructs momentum-resolved Hamiltonians $\( H(k_x, k_y) \)$ over a 2D Brillouin zone
- Computes eigenvalues and eigenvectors to extract band structures
- Calculates retarded Green’s functions $\( G(k, \omega) \)$, then integrates over k to obtain:
  - Local Green’s function $\( G_{\text{loc}}(\omega) \)$
  - Orbital- and cluster-resolved spectral functions and DOS
  - Total spectral weight $\( A(\omega) = -\frac{1}{\pi} \mathrm{Im} \, G_{\text{loc}} \)$
- Computes Matsubara Green’s functions $\( G(i\omega_n) \)$ from spectral representation
- Computes imaginary-time Green’s function $\( G(-it) \)$ via Fourier transform
- Tracks doping and local charge fluctuations as a function of modulation amplitude

Outputs:
- Spectral function plots and orbital-resolved DOS
- Cluster-resolved doping vs site index
- Doping data saved to: `doping_ilmu0.dat`
- Optional display of Green’s function in frequency and time domains

Requirements:
- numpy
- scipy
- matplotlib
- tqdm
- termcolor

Usage:
1. Edit physical parameters in the script (e.g., N, modulation amplitudes, hopping terms).
2. Run the script:
   python il16_y_modulation.py
3. Visualization and data files will be generated automatically.

Notes:
- The model assumes translational symmetry within N-mers and periodic boundary conditions.
- Code is computationally intensive; large values of Lx × Ly and fine frequency grids will impact performance.

