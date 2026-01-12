---
title: NEGF Implementation Quick Reference
---

# NEGF Implementation Quick Reference {#rt_negf_quick_ref}

**Concise summary of non-equilibrium Green's function capabilities in YAMBO real-time.**

---

## Direct Answers to Key Questions

### Q1: Are Occupations Evolved and Recomputed?

**YES - FULLY IMPLEMENTED**

$$f_n(k,t) = \text{Im}[G^<_{nn}(k,t,t)]$$

- Extracted from diagonal of Green's function at each timestep
- Evolved through Kadanoff-Baym equation
- Density $\rho(r,t)$ recomputed from occupations
- Hartree and XC potentials updated accordingly

Source: `src/real_time_propagation/RT_occupations_eval.F`

### Q2: Is Dielectric Function Updated in Real-Time?

**NO - USES STATIC EQUILIBRIUM SCREENING**

| Stage | Status | Details |
|-------|--------|---------|
| Compute $\varepsilon(q,\omega)$ | Once | Before RT simulation (GW) |
| Use in RT | Static | Fixed throughout propagation |
| Update with carriers | No | Screening does NOT adapt |

Mean-field self-energy:

$$\Delta\Sigma^{\text{xc}} = K_{\text{eq}} \cdot [G(t) - G_{\text{ref}}]$$

The kernel $K_{\text{eq}}$ is fixed. Dynamic $W(q,\omega,t)$ not computed.

Source: `src/collisions/COLLISIONS_compose_rt.F:86-87`

### Q3: Is BSE Kernel Updated in Real-Time?

**NO - BSE NOT IMPLEMENTED IN RT**

| Item | Status |
|------|--------|
| BSE kernel computation | Not called |
| Dynamic $K(q,t)$ | No time dependence |
| Excitonic renormalization | Missing |

The BSE module (`src/bse/`) exists but is never invoked during real-time propagation.

---

## Implementation Status Table

| Feature | Status | Type | Notes |
|---------|--------|------|-------|
| Occupation evolution | Implemented | Dynamic | Fully recomputed each step |
| Carrier creation (laser) | Implemented | Dynamic | Via external field coupling |
| Hartree potential | Implemented | Dynamic | Updated with carrier density |
| XC potential (TDDFT) | Implemented | Dynamic | Only in TDDFT mode |
| Collision kernel (TDSEX) | Implemented | Static | Equilibrium kernel, not updated |
| Band-filling effects | Implemented | Dynamic | Via occupation depletion |
| Nonlinear current | Implemented | Dynamic | From Green's function evolution |
| Scattering rates | Implemented | Dynamic | el-ph, el-el collision terms |
| | | | |
| Dielectric function $\varepsilon(q,\omega,t)$ | Missing | Static | Pre-computed at equilibrium |
| Plasma frequency shift | Missing | - | No carrier-dependent screening |
| Excitonic screening | Missing | - | No BSE kernel |
| BSE kernel $K(t)$ | Missing | - | Module not used in RT |
| Many-body renormalization | Missing | - | Static kernel approximation |
| Exciton-carrier interactions | Missing | - | Single-particle only |

---

## Screening Treatment Modes

### TDDFT Mode: Direct Functional Approach

$$H^{\text{RT}} = H^{\text{EQ}} + V^H[\rho(t)] + V^{\text{xc}}[\rho(t)] + H^{\text{ext}}(t)$$

**Properties**:
- Hartree potential: Dynamically updated
- XC potential: From functional (LDA/GGA)
- Screening: Implicit in XC functional
- Implementation: Full

**Advantages**:
- Hartree is truly dynamic
- Computationally efficient

**Limitations**:
- XC is semi-local (no non-local effects)
- Correlation screening incomplete

### TDSEX Mode: Static Kernel Approach

$$H^{\text{RT}} = H^{\text{EQ}} + K_{\text{eq}} \cdot [G(t) - G_{\text{ref}}] + H^{\text{ext}}(t)$$

**Properties**:
- Kernel from equilibrium GW
- Applied to density perturbation
- Screening: Static from GW

**Advantages**:
- Better band structure from GW
- Includes many-body effects

**Limitations**:
- Screening frozen in time
- Valid only for weak perturbations
- Kernel never recomputed

### TD-Hartree Mode: Classical Screening

$$H^{\text{RT}} = H^{\text{EQ}} + V^H[\rho(t)] + H^{\text{ext}}(t)$$

**Properties**:
- Hartree only (no exchange/correlation)
- Fully dynamic

**Advantages**:
- Fastest computation
- Exact for density waves

**Limitations**:
- No exchange effects
- No many-body correlations

---

## Timescale Validity

| Timescale | Carrier Dynamics | Screening Response | Status |
|-----------|-----------------|-------------------|--------|
| 1-100 fs | Fully captured | Static equilibrium | Valid |
| 100 fs - 1 ps | Captured with scattering | Static equilibrium | Valid for weak fields |
| Sub-fs (excitonic) | Not captured | Static equilibrium | Invalid |

The static screening approximation breaks down when screening changes significantly on the timescale of interest.

---

## Physical Effects: What is Included

| Effect | Captured | How |
|--------|----------|-----|
| Carrier creation by laser | Yes | External field coupling in H |
| Band filling (bleaching) | Yes | Occupation depletion |
| Coulomb repulsion | Yes | Hartree potential with carriers |
| Thermal effects | Partial | Via scattering rates |
| Momentum relaxation | Yes | Electron-phonon scattering |
| Dephasing | Yes | Collision self-energy |
| Nonlinear susceptibility | Yes | From Green's function evolution |

## Physical Effects: What is Missing

| Effect | Missing | Why |
|--------|---------|-----|
| Excitonic renormalization | Yes | BSE kernel not computed |
| Screening buildup | Yes | Dielectric function static |
| Plasma effects | Yes | Dynamic polarization not tracked |
| Many-body dephasing | Yes | Beyond collision integral |
| Band gap shift | Yes | Screening effects missing |
| Strong-field breakdown | Yes | Static kernel invalid at high fields |

---

## Code Organization Summary

### Dynamic (Updated Each Step)

- Occupations: `src/real_time_propagation/RT_occupations_eval.F`
- Density: `el_density_matrix()` called in `RT_Hamiltonian.F`
- Hartree potential: `V_Hartree()` in `RT_Hamiltonian.F:105`
- XC potential: `XC_potential_driver()` in `RT_Hamiltonian.F:119` (TDDFT mode)

### Static (Fixed During Propagation)

- Collision kernel: `src/collisions/COLLISIONS_compose_rt.F:55-95`
- Kernel application: Uses equilibrium `HXC_collisions` object
- Screening: Pre-computed W from GW setup

### Not Implemented

- Polarization function: `src/pol_function/` modules not called
- Dielectric response: No dynamic computation
- BSE kernel: `src/bse/` modules not invoked
- Time-dependent screening: Missing framework

---

## Recommended Practices

### For Accurate Results

1. **Choose appropriate mode**:
   - Weak fields: TDDFT-LDA or TD-Hartree
   - Better band structure: TDSEX (accept static screening)
   - Maximum simplicity: TD-Hartree

2. **Verify regime validity**:
   - Compare results at different field strengths
   - If results depend strongly on intensity, static kernel breaks down
   - Use convergence tests with different band ranges

3. **Physical interpretation**:
   - Carrier dynamics: Quantitatively accurate
   - Screening effects: Qualitative at best
   - Exciton physics: Not captured

### For Strong-Field Physics

1. Reduce field strength and extrapolate
2. Use TDDFT mode (more accurate than TDSEX)
3. Post-process with BSE at fixed carrier densities
4. Consider limitations in published results

---

## Comparison with Full NEGF

### Time-Diagonal Propagation (Implemented)

$$i\frac{\partial \rho(t)}{\partial t} = [H(t), \rho(t)]$$

- Propagate: $G^<(t,t)$ only (single time)
- Memory: $O(N_b^2 \cdot N_k)$ per time step
- Cost: $O(N_t)$ for $N_t$ time steps
- Available: All standard observables

### Full Two-Time Propagation (Not Implemented)

$$i\frac{\partial G^<(t,t')}{\partial t} = [H(t), G^<(t,t')]$$

- Would need: Two independent times $t$ and $t'$
- Memory: $O(N_t^2 \cdot N_b^2 \cdot N_k)$ - enormous
- Cost: $O(N_t^2)$ - prohibitively expensive
- Would enable: Memory effects, spectral functions

---

## References

**For Theory**: [NEGF Theory](YAMBO_RT_NEGF_THEORY.md)

**For Implementation Details**: [NEGF Implementation](YAMBO_RT_NEGF_IMPLEMENTATION.md)

**For Screening**: [Screening Theory](YAMBO_RT_SCREENING_THEORY.md)

**Source Code**:
- Main driver: `src/real_time_drivers/RT_driver.F`
- Hamiltonian: `src/real_time_hamiltonian/RT_Hamiltonian.F`
- Occupations: `src/real_time_propagation/RT_occupations_eval.F`
- Collisions: `src/collisions/COLLISIONS_compose_rt.F`

---

## Documentation Links

| Document | Purpose |
|----------|---------|
| [Structure Guide](\_DOCUMENTATION_STRUCTURE.md) | How to navigate documentation |
| [NEGF Theory](YAMBO_RT_NEGF_THEORY.md) | Theoretical foundations |
| [NEGF Implementation](YAMBO_RT_NEGF_IMPLEMENTATION.md) | What is actually coded |
| [Screening Theory](YAMBO_RT_SCREENING_THEORY.md) | Dielectric response details |
| [Screening Implementation](YAMBO_RT_SCREENING_IMPLEMENTATION.md) | How screening is computed |
