---
title: Exciton-Phonon Coupling Theory
---

# Exciton-Phonon Coupling Theory

## Overview

Exciton-phonon coupling describes the interaction between electronic excitations (excitons) and crystal vibrations (phonons). This coupling is essential for understanding several important physical phenomena:

- Non-radiative decay and exciton lifetimes
- Phonon-assisted optical transitions (photoluminescence and absorption)
- Excitonic linewidths and dephasing
- Temperature-dependent optical properties

The theoretical framework presented here follows the work by Paleari et al. and implements the exciton-phonon coupling matrix elements within the context of the Bethe-Salpeter equation for excitonic effects.

## Fundamental Concepts

### Exciton-Phonon Hamiltonian

The exciton-phonon interaction is described through a perturbation to the electronic Hamiltonian caused by atomic displacements. When atoms vibrate, the electronic potential changes, creating a coupling mechanism.

The electron-phonon matrix elements define how phonons scatter electrons:

$$
g_{n,m}(k,k+q,\nu,\alpha) = \langle \psi_n(k)|\nabla V_\alpha(q)|\psi_m(k+q)⟩ / \sqrt{(M_\alpha \omega_q,\nu)}
$$

where:
- $n, m$ are electronic band indices
- $k, q$ are electron and phonon wave vectors
- $\nu$ is the phonon branch index
- $\alpha$ labels the atomic species
- $\omega_q$, $\nu$ is the phonon frequency
- $M_\alpha$ is the atomic mass

### Exciton Wavefunctions

An exciton is a composite bound state of an electron in the conduction band and a hole in the valence band. The exciton wavefunction in the Tamm-Dancoff approximation is:

$$
|S,q\rangle = \sum_{k,v,c} A_{k}^{S}(q) c^\dagger_{c,k+q} c_{v,k} |GS⟩
$$

where:
- $S$ labels the excitonic state
- $k$ is the relative electron-hole momentum
- $c^\dagger$, $c$ are creation and annihilation operators
- $A_{k}^{S}(q)$ are expansion coefficients from the Bethe-Salpeter equation
- |GS⟩ is the ground state

## Exciton-Phonon Matrix Elements

### Definition

The exciton-phonon matrix element describes how a phonon couples two excitonic states. The **most general form** allows arbitrary exciton and phonon momenta:

$$
G^{\lambda, exc}_{\alpha,\beta}(Q,q) = ⟨\alpha,q+Q|\partial^\lambda_q V_|\beta,Q⟩
$$

where:
- $Q$ is the initial exciton wavevector 
- $q$ is the phonon wavevector
- $\alpha$ is the final exciton state at wavevector $Q+q$
- $\beta$ is the initial exciton state at wavevector $Q$
- $\nu$ is the phonon branch
- $\partial^\lambda_q V$ is the directional derivative of the total Kohn-Sham potential with respect to the phonon displacement vector of mode $\lambda$ and momentum $q$.

The code currently implements the **restricted case** where the initial exciton is fixed at the Gamma point ($Q=0$):

$$
G^{exc}_{\alpha,\beta}(q,\nu) = ⟨\alpha,q|\partial_q V_\nu|\beta,0⟩
$$


### Decomposition into Electron and Hole Contributions

The exciton-phonon matrix element can be decomposed into contributions from electron and hole scattering:

$$
G^{exc}_{\alpha,\beta}(q,ν) = \sum_{k,v,c,v',c'} A^*_{k+q,v',c'}(\beta) g_{c',c}(k,q,\nu) A_{k,v,c}(\alpha) + 
                    \sum_{k,v,c,v',c'} A^*_{k+q,v',c'}(\beta) g_{v,v'}(k,q,ν) A_{k,v,c}(\alpha)
$$

The first term represents phonon-induced electron scattering (electron contribution), and the second represents phonon-induced hole scattering (hole contribution).

The individual scattering matrix elements are:

$$
g_{c',c}(k,q,\nu) = ⟨c',k+q|\partial V_\nu|c,k⟩ / √(2\omega_q,\nu)
$$

$$
g_{v,v'}(k,q,\nu) = ⟨v,k|\partial V_\nu|v',k-q⟩ / √(2\omega_q,\nu)
$$
Note the opposite signs in the q-dependence: electrons pick up momentum from the phonon (k+q), while holes lose momentum (k-q).

## Phonon-Assisted Optical Transitions

### One-Phonon Processes

Phonon-assisted optical absorption allows transitions that would otherwise be forbidden by momentum conservation. The process can be decomposed as:

1. Virtual phonon absorption (or emission for PL)
2. Electronic transition between states with modified momentum
3. Recombination or absorption of the real photon

The contribution to the optical susceptibility from one-phonon processes:

$$
\chi_n = \sum_{\alpha,\beta,\nu} |G^{exc}_{\alpha,\beta}(q,\nu)|^2 / (E_\alpha - E_\beta - \omega_\nu - i\Gamma)
$$

The energy denominator shows that transitions occur near $E_\alpha - E_\beta = \omega_nu$ (phonon-assisted resonance).

### Satellite Peaks

In photoluminescence, one-phonon processes create satellite peaks at energies:

- $E_{sat}^{(-)} = E_{exc} - \omega_{ph}$  (phonon emission satellite)
- $E_{sat}^{(+)} = E_{exc} + \omega_{ph}$  (phonon absorption satellite)

The relative intensity depends on phonon occupations (Bose distribution) and exciton populations.

## Exciton Lifetime from Phonon Scattering

### Non-Radiative Lifetime

The non-radiative lifetime of an exciton due to phonon scattering is given by the imaginary part of the self-energy:

$$
1/\tau_{nr} = -2 Im[\Sigma(E)]
$$

The exciton self-energy from phonon scattering:

$$
\Sigma(E) = \sum_{\alpha,\nu} |G_{\alpha}(\nu)|^2 [n_\nu/(E - E_\alpha - \omega_\nu + i\Gamma) + (n_\nu+1)/(E - E_\alpha + \omega_\nu + i\Gamma)]
$$

where:
- $n_\nu$ is the Bose occupation of the phonon mode
- The first term represents phonon absorption (cooling process)
- The second term represents phonon emission (heating process)

### Temperature Dependence

The lifetime temperature dependence comes from the Bose distribution:

$$
n_\nu(T) = 1 / (exp(\omega_\nu/k_B T) - 1)
$$

At low temperature, phonon population approaches zero and non-radiative decay is suppressed. At high temperature, increased phonon population leads to faster non-radiative decay.

## Linewidth and Dephasing

### Homogeneous Linewidth

The exciton linewidth (FWHM) is related to the inverse lifetime:

$$
\Gamma = \hbar/\tau = -2 Im[\Sigma(E)]
$$

This describes the dephasing caused by interaction with the phonon bath. The Lorentzian lineshape:

$$
L(\omega) = \Gamma/(2\pi[(\omega - E)^2 + (\Gamma/2)^2])
$$

### Broadening Mechanisms

In the code, the linewidth includes contributions from:

1. Phonon scattering (main contribution computed here)
2. Other scattering mechanisms (radiative decay, defects) - treated as phenomenological damping

## Double-Grid Approach

For efficient computation, especially in systems with weak exciton binding, the code can employ a double-grid (coarse IBZ + refined fine grid) to:

1. Compute exciton-phonon matrix elements on a finer q-grid
2. Account for band structure variations across the Brillouin zone
3. Improve convergence for materials with small exciton binding energies


## References

The exciton-phonon implementation follows:

1. Paleari, F., "Phonon-assisted transitions in photoexcited semiconductors," PhD Thesis (2019)
2. Poncé, S., Margine, E., & Giustino, F., "Toward predictive many-body calculations of phonon-limited carrier mobilities in semiconductors," Physical Review B 97, 121201 (2018)
3. Standard theory as presented in Grosso & Parravicini "Solid State Physics" and Mahan "Many-Particle Physics"
