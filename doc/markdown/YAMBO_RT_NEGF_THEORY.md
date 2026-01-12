---
title: Non-Equilibrium Green's Functions in Real-Time Simulations
author: Based on YAMBO implementation analysis
date: 2024-12-04
---

# Non-Equilibrium Green's Functions - Theory {#rt_negf_theory}

**Objective**: Understand the theoretical foundations of non-equilibrium Green's functions as used in YAMBO real-time simulations.

**Prerequisites**: 
- [Kadanoff-Baym Equations](YAMBO_RT_KADANOFF_BAYM_THEORY.md)
- [Equilibrium Green's Functions](YAMBO_THEORY_DOCUMENTATION.md)

**Reading Time**: 20-30 minutes

---

## 1. Overview {#overview}

Non-equilibrium Green's functions (NEGF) provide a systematic framework for describing systems driven away from equilibrium by external perturbations. The approach is particularly valuable for:

- Ultrafast dynamics following laser excitation
- Carrier dynamics and thermalization
- Optical response and high-harmonic generation
- Time-resolved spectroscopy

The fundamental object is the **time-dependent Green's function**, which contains all information about the quantum state of the system.

---

## 2. Two-Time Green's Functions {#two_time_gf}

### Definition

The most general form of the lesser Green's function involves two independent time arguments:

$$G^<_{nm}(k,t,t') = i\langle c^\dagger_m(k,t') c_n(k,t) \rangle$$

where:
- $n, m$ are band indices
- $k$ is crystal momentum  
- $t, t'$ are two independent time variables
- $c^\dagger, c$ are creation and annihilation operators

### Physical Interpretation

The diagonal elements $G^<_{nn}(k,t,t) = if_n(k,t)$ represent occupations at time $t$. The off-diagonal elements $G^<_{nm}(k,t,t)$ with $n \neq m$ represent coherences between states.

For $t \neq t'$, the function carries information about:

- **Memory effects**: How the system remembers its past states
- **Correlations**: Time-time correlations in the density matrix  
- **Dephasing**: Decay of quantum coherences

### Two-Time Components

The complete Green's function has four independent components:

$$G^<(t,t') \quad G^>(t,t') \quad G^R(t,t') \quad G^A(t,t')$$

These satisfy the causality relation:

$$G^>(t,t') - G^<(t,t') = G^R(t,t') - G^A(t,t')$$

---

## 3. The Kadanoff-Baym Equations {#kbe}

### General Form

The complete set of Kadanoff-Baym equations for the two-time Green's function is:

$$i\frac{\partial G^<(t,t')}{\partial t} = [H(t), G^<(t,t')] + \int_{t'}^t d\tau\, \Sigma^R(t,\tau)G^<(\tau,t')$$

$$i\frac{\partial G^<(t,t')}{\partial t'} = -[H(t'), G^<(t,t')] - \int_{t}^{t'} d\tau\, G^<(t,\tau)\Sigma^A(\tau,t')$$

These are **integro-differential equations** with memory kernels in the self-energy.

### Physical Meaning

- First equation: Forward time evolution with history effects
- Second equation: Backward time evolution with history effects
- Memory kernels: Scattering processes with retardation

---

## 4. Time-Diagonal Approximation {#time_diagonal}

### Definition

The **time-diagonal approximation** reduces the two-time Green's function to a single-time density matrix:

$$G^<_{nm}(k,t,t) = \rho_{nm}(k,t)$$

This removes the $t'$ dependence and yields:

$$i\frac{\partial \rho(t)}{\partial t} = [H(t), \rho(t)] + \text{collision term}$$

This is a **single-time equation** that can be solved forward in time.

### Validity Conditions

The approximation is valid when:

1. **Short memory**: System does not depend on distant past
2. **Markovian scattering**: Collision self-energy is instantaneous
3. **Weak coupling**: Perturbation is small compared to system energy scale

### Implications

**What is preserved**:
- Occupations and populations
- Coherences at the same time: $\rho_{nm}(k,t)$
- Standard observables: current, energy, polarization

**What is lost**:
- Memory effects: $G^<(t,t')$ for $t \neq t'$
- Non-Markovian dynamics
- Spectral functions (require Fourier transform over $t-t'$)

---

## 5. Observables from NEGF {#observables}

### Carrier Populations

$$f_n(k,t) = \text{Im}[G^<_{nn}(k,t,t)] = \text{Im}[\rho_{nn}(k,t)]$$

The imaginary part of diagonal Green's function elements gives occupation factors.

### Electron Density

$$\rho(r,t) = -i\sum_{n,k} G^<_{nn}(k,t,t)|{\psi_{nk}(r)}|^2 = \sum_{n,k} f_n(k,t)|{\psi_{nk}(r)}|^2$$

Sum over all states weighted by occupations gives charge density.

### Current Density

$$J_\alpha(t) = -e\sum_{n,m,k} \text{Re}[\rho_{nm}(k,t)v^\alpha_{mn}(k)]$$

The current arises from both diagonal (intraband) and off-diagonal (interband) contributions.

### Energy

$$E_{\text{kin}}(t) = \sum_{n,k} f_n(k,t)E_n(k)$$

Kinetic energy from occupation-weighted eigenvalues.

---

## 6. Self-Energy Structure {#self_energy}

### General Decomposition

The total self-energy is split into coherent and collision parts:

$$\Sigma(t) = \Sigma^{\text{coh}}(t) + \Sigma^{\text{scatt}}(t)$$

### Coherent Self-Energy

$$\Sigma^{\text{coh}}_{\text{mean-field}} = \Delta \Sigma^{H}[G(t)] + \Delta \Sigma^{\text{xc}}[G(t)]$$

where:

- $\Delta \Sigma^H[G(t)]$ is the Hartree term (classical Coulomb repulsion)
- $\Delta \Sigma^{\text{xc}}[G(t)]$ is exchange-correlation effects

These can be computed as:

**Functional approach**:
$$\Delta \Sigma^{\text{xc}}[G(t)] = V^H[\rho(t)] + V^{\text{xc}}[\rho(t)]$$

Updated dynamically with the time-dependent density.

**Kernel approach**:
$$\Delta \Sigma^{\text{xc}} = K_{\text{eq}} \cdot [G(t) - G_{\text{ref}}]$$

where $K_{\text{eq}}$ is a static kernel computed at equilibrium.

### Collision Self-Energy

$$\Sigma^{\text{scatt}}(t) = \sum_\alpha \Sigma^\alpha_{\text{scatt}}(t)$$

Includes contributions from:

- Electron-electron scattering
- Electron-phonon coupling
- Dephasing processes
- Radiative recombination

---

## 7. Screening and Dielectric Response {#screening}

### Lindhard Polarization

The bare polarization function relates density fluctuations to perturbations:

$$\chi_0(q,\omega,t) = -\sum_{n,m,k} \frac{f_n(k,t) - f_m(k+q,t)}{E_m(k+q) - E_n(k) - \omega}$$

In equilibrium, this depends only on frequency. In real-time, it becomes time-dependent through the occupations $f_n(k,t)$.

### Dielectric Function

The full dielectric response is:

$$\varepsilon(q,\omega,t) = 1 - v(q)\chi(q,\omega,t)$$

where:

- $v(q) = 4\pi/(q^2)$ is the Coulomb kernel
- $\chi(q,\omega,t)$ includes both bubbles and vertex corrections

### Screened Coulomb Interaction

$$W(q,\omega,t) = \frac{v(q)}{\varepsilon(q,\omega,t)}$$

The screened interaction depends on the dynamic polarization response.

### Static vs. Dynamic Screening

**Static approximation**:
$$\varepsilon(q,t) = \varepsilon(q,\omega=0)$$

Use only the static response, independent of frequency and time.

**Dynamic approximation**:
$$\varepsilon(q,\omega,t) = 1 - v(q)\chi_0(q,\omega,t)$$

Full frequency and time dependence from polarization.

---

## 8. Mean-Field Approaches in NEGF {#mean_field}

### Independent Particle (IP)

$$H^{\text{RT}} = H^{\text{EQ}} + H^{\text{ext}}(t)$$

No mean-field corrections. Valid only for non-interacting systems.

### Time-Dependent Hartree (TD-Hartree)

$$H^{\text{RT}} = H^{\text{EQ}} + V^H[\rho(t)] + H^{\text{ext}}(t)$$

Includes classical Coulomb repulsion updated with carrier density. Exact for density oscillations.

### Time-Dependent DFT (TDDFT)

$$H^{\text{RT}} = H^{\text{EQ}} + V^H[\rho(t)] + V^{\text{xc}}[\rho(t)] + H^{\text{ext}}(t)$$

Hartree plus exchange-correlation effects. XC potential updated dynamically.

### Time-Dependent GW / Screened Exchange

$$H^{\text{RT}} = H^{\text{EQ}} + \Sigma^{\text{SEX}}[G_{\text{eq}}, W_{\text{eq}}] + H^{\text{ext}}(t)$$

Uses equilibrium screening. Provides better band structure but screening is static.

---

## 9. External Field Coupling {#external_field}

### Velocity Gauge

Coupling to electromagnetic field via vector potential:

$$H^{\text{ext}} = \frac{({\mathbf{p}} + {\mathbf{A}}(t))^2}{2m} - \frac{{\mathbf{p}}^2}{2m} = \frac{{\mathbf{A}}(t) \cdot {\mathbf{p}}}{m} + \frac{{\mathbf{A}}^2(t)}{2m}$$

The matrix elements are:

$$H^{\text{ext}}_{nm} = \frac{{\mathbf{A}}(t) \cdot {\mathbf{p}}_{nm}}{m} - \delta_{nm}\frac{{\mathbf{A}}^2(t)}{2m}$$

where ${\mathbf{p}}_{nm} = \langle n|{\mathbf{p}}|m\rangle$ is the momentum matrix element.

### Length Gauge

In the dipole approximation, the coupling is:

$$H^{\text{ext}} = {\mathbf{E}}(t) \cdot {\mathbf{r}}$$

Matrix elements:

$$H^{\text{ext}}_{nm} = {\mathbf{E}}(t) \cdot {\mathbf{r}}_{nm}$$

where ${\mathbf{r}}_{nm} = \langle n|{\mathbf{r}}|m\rangle$ is the position matrix element.

### Gauge Equivalence

The two formulations are physically equivalent but lead to different numerical behavior and convergence properties.

---

## 10. Applications and Limitations {#applications}

### Suitable Applications

1. **Carrier dynamics**: Thermalization, recombination timescales
2. **Nonlinear optical response**: Harmonic generation, high-field effects
3. **Transient spectroscopy**: Pump-probe measurements  
4. **Coherent control**: Pulse shaping optimization

### Current Limitations

1. **Markovian approximation**: No memory in scattering during propagation
2. **Static screening**: Screening kernels not updated in real-time
3. **Time-diagonal only**: Full two-time correlations not propagated
4. **Single-particle description**: Excitonic effects approximate

### Beyond NEGF Time-Diagonal

For phenomena requiring full two-time Green's functions:

- Non-Markovian dynamics
- Memory effects in scattering
- Spectral functions and spectroscopy
- Strong-field non-perturbative regimes

---

## See Also

- [Kadanoff-Baym Implementation](YAMBO_RT_KADANOFF_BAYM_IMPLEMENTATION.md): How KBE are solved
- [Occupation Evolution](YAMBO_RT_OCCUPATION_EVOLUTION_THEORY.md): Carrier dynamics theory
- [Screening Theory](YAMBO_RT_SCREENING_THEORY.md): Dielectric response details
- [NEGF Implementation](YAMBO_RT_NEGF_IMPLEMENTATION.md): What is actually implemented

## References

1. H. Haug and A.-P. Jauho, *Quantum Kinetics in Transport and Optics of Semiconductors* (Springer, 2008)
2. G. Stefanucci and R. van Leeuwen, *Nonequilibrium Many-Body Theory of Quantum Systems* (Cambridge, 2013)
3. L.V. Keldysh, *Sov. Phys. JETP* **20**, 1018 (1965)
4. L.P. Kadanoff and G. Baym, *Quantum Statistical Mechanics* (Benjamin, 1962)
5. C. Attaccalite et al., *Phys. Rev. B* **84**, 245110 (2011) - YAMBO real-time implementation
