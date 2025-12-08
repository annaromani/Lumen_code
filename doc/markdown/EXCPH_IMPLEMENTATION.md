---
title: Exciton-Phonon Coupling Implementation
---

# Exciton-Phonon Coupling Implementation

## Overview

This document describes how exciton-phonon coupling is computed in the LUMEN code. The implementation consists of several stages: coupling matrix element evaluation, optical property calculation, and lifetime computation.

## Workflow

The exciton-phonon calculation proceeds in the following order:

1. Load electron-phonon (gkkp) matrix elements from pre-computed databases
2. Rotate exciton wavefunctions to different q-points in the Brillouin zone
3. Compute exciton-phonon matrix elements by combining electron-phonon and exciton data
4. Calculate optical properties including phonon-assisted absorption
5. Compute exciton lifetimes and linewidths

## Computing Exciton-Phonon Matrix Elements

### Entry Point: EXCPH_gkkp_driver

The calculation is orchestrated by `EXCPH_gkkp_driver`. This routine:

1. Loads pre-computed electron-phonon matrix elements (gkkp) from disk
2. Allocates data structures for storing exciton-phonon results
3. Performs a q-point loop where for each q:
   - Exciton wavefunctions are rotated to the new q-point using symmetries
   - Individual phonon branches are processed
   - Results are written to disk

### Key Algorithm: EXCPH_gkkp_eval

The core calculation is performed by `EXCPH_gkkp_eval`, which implements the formula:

$$
G^{exc}_{\alpha,\beta}(q,\nu) = \sum_H A^*_{H}(\beta) [\xi_H^{hole}(q,\nu) - \xi_H^{elec}(q,\nu)] A_{H}(\alpha)
$$

The calculation proceeds in two steps:

#### Step 1: Compute $\xi$ Terms

For each exciton state $\beta$ (final state), compute the particle-phonon coupling:

$$
ξ_H^{hole}(q,\nu) = \sum_{v'} g_{c,v'}(k,q,\nu) A^*_{H'}(k,v',c,\beta)
$$

$$
ξ_H^{elec}(q,\nu) = \sum_{c'} g_{c',c}(k+q,q,ν) A^*_{H'}(k+q,v,c',\beta)
$$

where:
- H indexes an electron-hole pair (k, v, c)
- g terms are the electron-phonon matrix elements normalized by $\sqrt(2\omega_\nu)$

The code handles two contributions:
- Hole contribution (DOWN): electron stays fixed, hole scatters to v'
- Electron contribution (UP): hole stays fixed, electron scatters to c'

The minus sign reflects opposite momentum transfers for electrons vs holes.

#### Step 2: Compute Gkkp

For each exciton state α (initial state), sum over all electron-hole pairs:

$$
G^{exc}_{\alpha,\beta}(q,\nu) = \sum_H A_{H}(\alpha) [\xi_H^{hole} - \xi_H^{elec}]
$$

This involves conjugating the exciton wavefunction coefficient for state β and summing the Xi terms weighted by state α coefficients.

### Wavefunction Rotation

Exciton wavefunctions must be rotated from the reference $q=0$ to other q-points using:

$$
A_{k,v,c}(\alpha,q) = U^\dagger_{v,v'}(k_{sym},q_{sym}) U_{c,c'}(k+q_{sym},q_{sym}) A_{k,v,c}(\alpha)
$$

The unitary matrices come from diagonalizing the electronic Hamiltonian at each q-point and accounting for crystal symmetries via phase factors stored in ELEC_WF_phases.

## Computing Optical Properties

### Phonon-Assisted Absorption: EXCPH_X_phass

This routine calculates the optical susceptibility including phonon-assisted processes:

#### Algorithm

For each frequency $\omega$ and q-point:

1. Compute absorption resonance from direct exciton transitions
2. Add one-phonon satellite contributions from phonon absorption/emission
3. Include two-phonon contributions if necessary

The susceptibility includes three components:

$$
\chi_\alpha(\omega) = \chi_{phonon}(\omega) + \chi_{no-satellite}(\omega) + \chi_{satellite}(\omega)
$$

where:
- $\chi_{phonon}$: Indirect transitions where a photon combines with phonon creation/annihilation
- $\chi_{no-satellite}$: Direct exciton transitions with phonon broadening
- $\chi_{satellite}$: Pure phonon-assisted peaks at $E \pm \omega_\nu$

#### Key Formula

The phonon-assisted absorption contribution near frequency ω:

$$
\chi_\alpha^{ph-ass}(\omega) = \sum_{\beta,\nu} |G^{exc}_{\alpha,\beta}(q,\nu)|^2 / [(E_\beta - \omega) - \omega_\nu]
$$

This shows the resonance occurs at $\omega = E_\beta - \omega_\nu$ (phonon emission satellite).

### Luminescence (PL): EXCPH_X_phass

Photoluminescence involves exciton recombination with phonon assistance. The code computes:

1. Direct PL: Exciton radiates photon directly
2. Satellite PL: Exciton emits phonon then radiates
3. Temperature-dependent populations using Bose and Fermi distributions

The output distinguishes three contributions:
- Im[eps_ph]: Contribution from phonon-assisted processes
- Im[eps_nos]: No-satellite direct contribution  
- Im[eps_sat]: Pure satellite (resonant phonon-assisted) contribution

## Computing Exciton Lifetimes

### EXCPH_lifetime Algorithm

The non-radiative lifetime is computed from the exciton self-energy:

$$
1/\tau_{nr} = 2 \sum_{\beta,\nu} |G^{exc}_{\alpha,\beta}(q,\nu)|^2 \delta(E_\alpha - E_\beta \mp \omega_\nu)
$$

The delta function is replaced by a Lorentzian with broadening parameter `L_damping`.

#### Implementation Details

For each initial exciton state α:

1. Loop over all final states β
2. For each phonon branch ν
3. Check if energy conservation is approximately satisfied: $|E_\alpha - E_\beta - \omega_\nu| < \textrm{broadening}$
4. Accumulate the lifetime contribution: $|G^{exc}|^2 × \textrm{density-of-states}$ factor

Two processes contribute:
- Phonon absorption ($\beta \rightarrow \alpha - \omega_\nu$): Creates phonon
- Phonon emission ($\beta → \alpha + \omega_\nu$): Annihilates phonon

Phonon populations at finite temperature are included:

$$
1/\tau_{nr} = \sum_{\beta,\nu} |G^{exc}_{\alpha,\beta}|^2 [n_\nu/(E_\alpha-E_\beta-\omega_\nu) + (n_\nu+1)/(E_\alpha-E_\beta+\omega_\nu)]
$$

## Output Data Structures

### EXCPH_gkkp Array

Stores the computed exciton-phonon matrix elements:

EXCPH_gkkp(il, i_beta, i_alpha)

Dimensions:
- il: Phonon branch (1 to n_branches)
- i_beta: Final exciton state
- i_alpha: Initial exciton state

### EXCPH_gkkp_sq Array

Squared magnitude of coupling elements:

EXCPH_gkkp_sq(il, i_beta, i_alpha) = |EXCPH_gkkp(il, i_beta, i_alpha)|^2

Used for lifetime calculations and diagnostics.

### EXC_LifeTime Array

Non-radiative lifetime in atomic units (Hartree^{-1}):

EXC_LifeTime(i_alpha)

Inverse gives the decay rate in energy units.

## Computational Efficiency

### Double-Grid Implementation

For materials with weak binding, exciton-phonon coupling is computed on a finer q-grid (double-grid) to:

1. Resolve fine structure in the phonon coupling
2. Improve convergence for small-gap materials
3. Account for band structure curvature effects

The double-grid setup uses interpolation to map between coarse and fine grids.

### Parallelization

The calculation distributes work over MPI processes:

- Q-point parallelization: Different q-points computed on different ranks
- Phonon branch parallelization: Different branches handled independently
- Band-pair parallelization: Loop over electron-hole pairs distributed

### Memory Optimization

Several strategies reduce memory footprint:

1. Sequential I/O of electron-phonon matrix elements (not stored entirely in memory)
2. On-the-fly computation of unitary transformations for wavefunction rotation
3. Temporary storage of intermediate results with immediate deletion

## Options and Flags

### Configuration Parameters

Key input parameters controlling the calculation:

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| ConstElph | Logical | False | Use constant (q,ν-independent) coupling |
| AbsElph | Logical | False | Use absolute value of coupling |
| AbsExc | Logical | False | Use absolute value of wavefunctions |
| ElectronContributionOnly | Logical | False | Include only electron scattering |
| HoleContributionOnly | Logical | False | Include only hole scattering |
| KeepDeg | Logical | False | Maintain degenerate state splitting |
| NoMatrxEl | Logical | False | Skip matrix element computation |
| DbGdWEIGHTs | Logical | False | Use weights on double-grid |
| DbGdOnlyPh | Logical | False | Apply double-grid only to phonons |
| SatFromHPert | Logical | False | Use perturbation theory for satellites |
| ExcPhOffDiago | Logical | False | Include off-diagonal lifetime terms |

## Validation and Debugging

### Consistency Checks

The code includes checks for:

1. Correspondence between input and output exciton states
2. Energy conservation near resonances
3. Wavefunction normalization
4. Symmetry operations consistency

### Output Files

Generates several output files for analysis:

- excph_gkkp: Exciton-phonon matrix elements database
- lifetimes: Non-radiative lifetimes per exciton state
- eps_bse_ph_ass: Phonon-assisted dielectric function
- pl_bse_ph_ass: Phonon-assisted luminescence spectrum

## Performance Considerations

### Computational Cost

The main computation scales as:

O(N_q × N_phonon × N_exc^2 × N_kpoints × N_bands^2)

where:
- N_q: Number of q-points
- N_phonon: Number of phonon branches
- N_exc: Number of exciton states
- N_kpoints: Number of k-points in band structure
- N_bands: Number of bands

For typical systems (100-1000 q-points, 3 phonon branches, 10-100 excitons), computation takes minutes to hours on modern hardware.

## References

Implementation details follow:

1. Paleari, F., et al., "Exciton-phonon coupling in bulk crystals," Physical Review B (2021)
2. Salpeter, E.A., "Mass corrections to the fine structure of hydrogen-like atoms," Physical Review 87, 328 (1952)
3. Code comments in EXCPH_gkkp_eval.F and related subroutines
