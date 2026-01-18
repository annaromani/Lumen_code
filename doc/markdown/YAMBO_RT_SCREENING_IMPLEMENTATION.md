---
title: Screening and Dielectric Response - Implementation
---

# Screening Implementation {#rt_screening_implementation}

**Purpose**: Describe how screening is computed and applied in YAMBO real-time simulations.

**Source Code**: Collision kernels, GW setup, real-time Hamiltonian

**Complexity**: Intermediate

---

## 1. Overview {#overview}

YAMBO implements screening through two complementary strategies:

1. **Hartree screening**: Computed dynamically at each real-time step from carrier density
2. **GW/BSE screening**: Pre-computed at equilibrium and applied statically during propagation

This document describes both approaches, their implementations, and limitations.

---

## 2. Hartree Screening - Dynamic Implementation {#hartree_screening}

### Hartree Potential Computation

**File**: `src/real_time_hamiltonian/RT_Hamiltonian.F:99-107`

The Hartree potential is computed from the charge density:

```fortran
allocate(delta_rho(fft_size))
call el_density_matrix(dG_in, E, k, delta_rho, RT_bands(1))
rho_n = rho_reference + delta_rho

if (eval_HARTREE) call V_Hartree(delta_rho, V_hartree_sc)
```

The procedure:

1. **Extract density**: Compute $\rho(r,t) = \rho_{\text{ref}} + \Delta\rho(t)$ from occupations
2. **Solve Poisson**: Compute $V^H(r,t)$ from total density
3. **Store potential**: Save for matrix element computation

### Poisson Solver

**Subroutine**: `V_Hartree()`

Solves Poisson equation in reciprocal space:

$$V^H(q,t) = \frac{4\pi}{q^2} \rho(q,t)$$

Implementation uses:
- FFT to transform density to Fourier space
- Division by $|q|^2$
- Inverse FFT back to real space

### Matrix Elements

**File**: `src/real_time_hamiltonian/RT_Hamiltonian.F:151`

Hartree self-energy matrix elements:

```fortran
call V_real_space_to_H(ik,i_sp_pol,H_nl_sc,WF,'def+xc',&
                       V=V_hartree_sc,Vxc=V_xc_sc)
```

Computes:

$$[\Delta\Sigma^H]_{nm}(k,t) = \langle \psi_n|V^H(r,t)|\psi_m \rangle$$

This is **fully dynamic**: updated at each time step with current density.

### Hartree Density Response

The Hartree screening provides a classical Coulomb response. For a uniform electron gas:

$$\varepsilon_H(q,t) \approx 1 + \frac{\omega_p^2(t)}{q^2 v_s^2}$$

with $\omega_p(t) \propto \sqrt{n(t)}$.

---

## 3. XC Screening - Conditional Implementation {#xc_screening}

### TDDFT Mode: Dynamic XC

**File**: `src/real_time_hamiltonian/RT_Hamiltonian.F:118-121`

```fortran
if (eval_DFT.and..not.l_use_Hxc_collisions) then
  call XC_potential_driver(E,k,WF_kind,WF_xc_functional,1)
  V_xc_sc = V_xc - REF_V_xc_sc
endif
```

The XC potential is computed from current density at each step.

**Implementation**: Fully dynamic, updated each time step.

### Functionals Supported

| Functional | Type | Screening | Status |
|-----------|------|-----------|--------|
| LDA | Local | Implicit | Implemented |
| GGA | Semi-local | Implicit | Implemented |

---

## 4. Collision Kernel - Static Implementation {#collision_kernel}

### Kernel Application in RT

**File**: `src/collisions/COLLISIONS_compose_rt.F:55-95`

At each time step, the static kernel is applied:

```fortran
RT_Vnl_xc(i_n,i_m,i_k,i_sp) = RT_Vnl_xc(...) &
  - cI * HXC_COLL_element(...)%v_c(i_comp) &
  * dG_avg(ib,ibp,i_kmq,ip_sp)
```

The self-energy is:

$$\Sigma^{\text{xc}}_{\text{nm}}(k,t) = \sum K_{nm,ij}^{\text{eq}} \cdot G^<_{ij}(t)$$

where $K^{\text{eq}}$ is the **fixed equilibrium kernel** that is never updated during real-time propagation.

### Why Kernel is Static

The kernel is computed once at equilibrium and never updated because:

1. Computational cost would be prohibitive
2. Approximation valid for weak excitations
3. Implementation complexity

---

## 5. Screening Modes Comparison {#modes}

### Available Modes

| Mode | Formula | Screening | Implementation |
|------|---------|-----------|-----------------|
| TD-Hartree | $H^{\text{EQ}} + V^H[\rho(t)] + H^{\text{ext}}$ | Dynamic Hartree | Full |
| TDDFT-LDA | $H^{\text{EQ}} + V^H[\rho(t)] + V^{\text{xc}}[\rho(t)] + H^{\text{ext}}$ | Dynamic H + XC | Full |
| TDSEX | $H^{\text{EQ}} + K^{\text{eq}}[G(t)-G_{\text{ref}}] + H^{\text{ext}}$ | Static GW | Kernel fixed |

### Performance Comparison

| Aspect | TD-Hartree | TDDFT | TDSEX |
|--------|-----------|-------|-------|
| Speed | Fastest | Fast | Moderate |
| Hartree | Dynamic | Dynamic | Static |
| XC/Exchange | None | Dynamic | Static |
| Valid regime | Weak all fields | Weak fields | Very weak fields |

---

## 6. What is NOT Implemented {#missing}

### No Dynamic Dielectric Function

$\varepsilon(q,\omega,t)$ is never recomputed during real-time propagation.

### No Plasma Frequency Tracking

Plasma frequency shift $\omega_p(t) \propto \sqrt{n(t)}$ is not explicitly computed.

### No Excitonic Renormalization

Binding energy changes from screening are not included.

---

## See Also

- [Screening Theory](YAMBO_RT_SCREENING_THEORY.md): Physical background
- [NEGF Implementation](YAMBO_RT_NEGF_IMPLEMENTATION.md): Overall RT implementation
- [NEGF Quick Reference](YAMBO_RT_NEGF_QUICK_REFERENCE.md): Implementation status summary

## References

1. C. Attaccalite et al., *Phys. Rev. B* **84**, 245110 (2011)
2. Source code: `src/real_time_hamiltonian/RT_Hamiltonian.F`, `src/collisions/COLLISIONS_compose_rt.F`
