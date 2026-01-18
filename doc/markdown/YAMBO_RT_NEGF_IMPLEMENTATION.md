---
title: Non-Equilibrium Green's Functions - Implementation Status
---

# NEGF Implementation Analysis {#rt_negf_implementation}

**Purpose**: Assess what non-equilibrium Green's function features are actually implemented in YAMBO real-time simulations.

**Source Code**: Real-time drivers, Hamiltonian construction, collision kernels

**Complexity**: Intermediate

---

## 1. Overview {#overview}

This document provides a detailed analysis of NEGF implementation in yambo_rt, comparing theoretical ideals with actual code capabilities. The analysis is based on direct code review of the YAMBO source.

### Summary Table

| Feature | Implemented | Status | Evidence |
|---------|-------------|--------|----------|
| Two-time KBE | Time-diagonal only | Partial | `src/real_time_drivers/RT_driver.F` |
| Occupation evolution | Yes | Full | `RT_occupations_eval.F` |
| Occupation recomputation | Yes | Full | Main RT loop |
| Hartree potential update | Yes | Full | `RT_Hamiltonian.F` |
| XC potential update | Conditional | Partial | TDDFT mode only |
| Dynamic screening | No | Missing | Static kernel used |
| Dielectric function update | No | Missing | Pre-computed at equilibrium |
| BSE kernel evolution | No | Missing | BSE module not called |

---

## 2. Occupation Evolution - Full Implementation {#occupation_evolution}

### Extraction from Green's Function

**File**: `src/real_time_propagation/RT_occupations_eval.F:52-59`

Occupations are extracted from the diagonal elements of the time-dependent Green's function:

```fortran
RT_el_occ%N(ib,ik,i_sp) = aimag(G_rot(ib,ib))
RT_ho_occ%N(ib,ik,i_sp) = -aimag(G_rot(ib,ib)) + spin_occ
RT_el_occ%dN(ib,ik,i_sp) = aimag(dG(ib,ib,ik,i_sp))
```

This implements:

$$f_n(k,t) = \text{Im}[G^<_{nn}(k,t,t)]$$

where $G^<_{nn}$ is the diagonal lesser Green's function.

### Recalculation Each Time Step

**Location**: `src/real_time_drivers/RT_driver.F` - main time loop

At every time step:

1. Build Hamiltonian: `call RT_Hamiltonian(dG_in, A_input, E, k)`
2. Integrate: `call RT_Integrator(dG_in, E, k, q)`
3. Extract occupations: `call RT_occupations_eval(dG)`
4. Compute observables: `call RT_Observables(E, k, OBSERVABLES)`

The occupations are therefore recomputed fresh from the evolved Green's function at each step. This is fully implemented.

### Evolution Mechanism

Occupations evolve through three mechanisms:

**1. Laser-driven transitions** via external field coupling in Hamiltonian:

$$i\frac{\partial f_n(k,t)}{\partial t}\bigg|_{\text{laser}} \propto \text{Im}[H^{\text{ext}}_{nm} \rho_{mn}(t)]$$

**2. Hartree repulsion** from carrier density:

$$i\frac{\partial f_n(k,t)}{\partial t}\bigg|_H = \text{Im}[V^H[\rho(t)] \cdot G^<(t)]$$

**3. Scattering** from collision self-energy:

$$i\frac{\partial f_n(k,t)}{\partial t}\bigg|_{\text{scatt}} = -[\Sigma^{\text{scatt}}(t)]_{nn}$$

All three mechanisms are implemented.

---

## 3. Hartree Potential Update - Full Implementation {#hartree_update}

### Density Computation

**File**: `src/real_time_hamiltonian/RT_Hamiltonian.F:99-101`

```fortran
allocate(delta_rho(fft_size))
call el_density_matrix(dG_in, E, k, delta_rho, RT_bands(1))
rho_n = rho_reference + delta_rho
```

The electron density is computed from occupations:

$$\rho(r,t) = \sum_{n,k} f_n(k,t)|{\psi_{nk}(r)}|^2$$

This includes all k-points and bands in the active space.

### Hartree Potential Computation

**File**: `src/real_time_hamiltonian/RT_Hamiltonian.F:105`

```fortran
if (eval_HARTREE) call V_Hartree(delta_rho, V_hartree_sc)
```

Solves Poisson equation:

$$V^H(r,t) = \int \frac{\rho(r',t)}{|{\mathbf{r}} - {\mathbf{r}}'|} d^3r'$$

The Hartree potential is:

- Computed from current density at each step
- Fully dynamic - adapts to carrier population changes
- Included in the RT Hamiltonian via matrix elements

### Matrix Elements

**File**: `src/real_time_hamiltonian/RT_Hamiltonian.F:151`

```fortran
call V_real_space_to_H(ik,i_sp_pol,H_nl_sc,WF,'def+xc',&
                       V=V_hartree_sc,Vxc=V_xc_sc)
```

Real-space potential is transformed to band basis:

$$[\Delta\Sigma^H]_{nm} = \langle \psi_n|V^H[\rho(t)]|\psi_m\rangle$$

This is fully implemented and updated at each time step.

---

## 4. XC Potential Update - Conditional Implementation {#xc_update}

### When XC is Updated: TDDFT Mode

**File**: `src/real_time_hamiltonian/RT_Hamiltonian.F:118-121`

```fortran
if (eval_DFT.and..not.l_use_Hxc_collisions) then
  call XC_potential_driver(E,k,WF_kind,WF_xc_functional,1)
  V_xc_sc = V_xc - REF_V_xc_sc
endif
```

When not using collision kernels, the XC potential is:

- Recomputed at each time step from current density
- Obtained from XC functional (LDA/GGA)
- Included in Hamiltonian as matrix elements

Implementation: **FULL** for TDDFT mode

$$\Delta\Sigma^{\text{xc}}(t) = V^{\text{xc}}[\rho(t)] - V^{\text{xc}}_{\text{ref}}[\rho_{\text{ref}}]$$

### When XC is NOT Updated: Collision Mode

**File**: `src/real_time_hamiltonian/RT_Hamiltonian.F:170-191`

```fortran
if (l_use_Hxc_collisions) then
  call Build_dG_avg(dG_in, dG_avg)
  call COLLISIONS_compose_rt(dG_avg)
  
  do i_sp_pol=1,n_sp_pol
    do ik=1,nXkibz
      if (l_H_contains_HARTREE) then
        Ho_plus_Sigma(:,:,ik_RT,i_sp_pol) = Ho_plus_Sigma(...) + RT_Vnl_xc(:,:,ik,i_sp_pol)
      else
        Ho_plus_Sigma(:,:,ik_RT,i_sp_pol) = RT_Vnl_xc(:,:,ik,i_sp_pol)
      endif
    enddo
  enddo
endif
```

When using collision kernels (TDSEX/TDCOHSEX mode), the XC self-energy is:

$$\Delta\Sigma^{\text{xc}}(t) = K_{\text{eq}} \cdot [G(t) - G_{\text{ref}}]$$

where $K_{\text{eq}}$ is the **equilibrium kernel** from GW calculation.

Implementation: **KERNEL IS STATIC** - not recomputed

---

## 5. Screening Kernel - Static Implementation {#screening_kernel}

### How Collision Kernels Are Applied

**File**: `src/collisions/COLLISIONS_compose_rt.F:55-95`

```fortran
do i_coll=1,HXC_collisions%N
  
  i_k = HXC_collisions%state(i_coll,3)
  i_n = HXC_collisions%state(i_coll,1)
  i_m = HXC_collisions%state(i_coll,2)
  
  RT_Vnl_xc(i_n,i_m,i_k,i_sp) = RT_Vnl_xc(...) &
    - cI * HXC_COLL_element(i_coll_mem)%v_c(i_comp) &
    * dG_avg(ib,ibp,i_kmq,ip_sp)
    
enddo

call PP_redux_wait(RT_Vnl_xc, COMM=...)
```

The self-energy is built as:

$$\Sigma^{\text{xc}}_{\text{nm}}(k,t) = -i \sum_{\text{collision}} K_{\text{nm}} \cdot G^<(t)$$

where $K_{\text{nm}}$ is the equilibrium collision kernel.

### Why the Kernel is Static

The collision kernel $K_{\text{eq}}$ is:

1. Computed once at equilibrium using GW
2. Stored in database: `HXC_collisions` object
3. Applied to every time step with the same values
4. Never updated based on time-dependent density or Green's function

This is an **approximation valid only for small perturbations** where $[G(t) - G_{\text{ref}}] \ll G_{\text{ref}}$.

### Averaging Procedure

**File**: `src/collisions/COLLISIONS_compose_rt.F:128-175`

```fortran
subroutine Build_dG_avg(dG_in, dG_avg)
  do i_k=1,nXkibz
    if (.not.PAR_IND_Xk_ibz%element_1D(i_k)) cycle
    
    do i_k_RT=RTibz%k_range(i_k,1), RTibz%k_range(i_k,2)
      dG_avg(:,:,i_k,i_sp) = dG_avg(:,:,i_k,i_sp) + dG_in(:,:,i_k_RT,i_sp)
    enddo
    
    dG_avg(:,:,i_k,i_sp) = dG_avg(:,:,i_k,i_sp) / real(RTibz%N_mapped(i_k),SP) / real(PAR_COM_Q_INDEX%n_CPU)
    
  enddo
end subroutine
```

The Green's function is averaged over equivalent k-points before kernel application.

---

## 6. Dielectric Function - Not Updated {#dielectric_function}

### Status: NOT IMPLEMENTED

The dielectric function $\varepsilon(q,\omega,t)$ is **never computed during real-time propagation**.

### Evidence

**Polarization code location**: `src/pol_function/`

Directories exist for dielectric response computation but are:
- Only called during GW initialization
- Not invoked in RT driver
- Not part of RT time loop

**Verification**: Searching RT drivers for dielectric function calls returns zero results.

### What Would Be Needed

To implement dynamic screening:

1. **Compute polarization at each step**:
   
$$\chi(q,\omega,t) = -\sum_{nm} \frac{f_n(k,t) - f_m(k+q,t)}{E_m(k+q) - E_n(k) - \omega}$$

2. **Update dielectric function**:
   
$$\varepsilon(q,\omega,t) = 1 - v(q)\chi(q,\omega,t)$$

3. **Recompute screened interaction**:
   
$$W(q,\omega,t) = \frac{v(q)}{\varepsilon(q,\omega,t)}$$

4. **Rebuild collision kernels**:
   
$$K(q,t) = \text{f}(W(q,t))$$

This would be **very expensive** (FFT-like scaling in frequency) and is not implemented.

---

## 7. BSE Kernel Evolution - Not Implemented {#bse_kernel}

### Status: NO BSE IN REAL-TIME

The Bethe-Salpeter equation kernel $K(t)$ is **never computed during real-time propagation**.

### Evidence

**BSE module location**: `src/bse/`

Files in this directory:
- `K_kernel.F` - BSE kernel
- `K_correlation_kernel_std.F` - Correlation part
- `K_screened_interaction.F` - Screening part
- Many others (20+ files)

**None of these are called** in the RT driver (`src/real_time_drivers/RT_driver.F`).

### BSE in Equilibrium

In equilibrium calculations:

1. Compute screened interaction $W(q,\omega=0)$
2. Build BSE kernel: $K = K_{\text{direct}} + K_{\text{exchange}}$
3. Diagonalize to get exciton energies and wavefunctions

### Missing: Dynamic BSE

For real-time exciton dynamics, would need:

1. **Compute W(q,t)** at each time step from dynamic screening
2. **Rebuild K(t)** from dynamic screening and vertex corrections
3. **Solve coupled e-h dynamics** including excitonic effects
4. **Track binding energy renormalization** from carrier screening

This is **not implemented** in yambo_rt.

### Closest Implementation

The only many-body correlation effects come from:

- Hartree mean field (classical screening)
- XC functional (if using TDDFT)
- Collision kernel (if using TDSEX)

All of these use **equilibrium or static** kernels, not dynamic ones.

---

## 8. Comparison: Theory vs. Implementation {#comparison}

### What Theory Predicts

**Full NEGF with dynamic screening**:

$$\Sigma^{\text{xc}}(t) = W(q,t) \cdot G(q,t)$$

where $W$ depends on the dynamic polarization from current carriers.

### What Is Actually Computed

**Static kernel approximation used in TDSEX**:

$$\Sigma^{\text{xc}}(t) = K_{\text{eq}} \cdot [G(t) - G_{\text{ref}}]$$

where $K_{\text{eq}} = W_{\text{eq}} - W_{\text{eq,bare}}$ is frozen at equilibrium.

### Validity Regime

The static kernel approximation is valid when:

$$|G(t) - G_{\text{ref}}| \ll |G_{\text{ref}}|$$

Typically requires:

- Weak excitation (small optical fields)
- Short interaction times (before significant excitation)
- Low carrier densities (less screening change)

---

## 9. Code Locations Reference {#code_locations}

### Dynamic Components

| Feature | File | Lines | Status |
|---------|------|-------|--------|
| Occupation extraction | `RT_occupations_eval.F` | 30-75 | Full |
| Density computation | `el_density_matrix()` | Various | Full |
| Hartree potential | `V_Hartree()` | Various | Full |
| Hartree in Hamiltonian | `RT_Hamiltonian.F` | 105-151 | Full |
| XC potential (TDDFT) | `XC_potential_driver()` | Various | Full |

### Static Components

| Feature | File | Lines | Status |
|---------|------|-------|--------|
| Collision kernel | `COLLISIONS_compose_rt.F` | 55-95 | Static |
| Kernel application | `COLLISIONS_compose_rt.F` | 86-87 | Static |
| dG averaging | `COLLISIONS_compose_rt.F` | 128-175 | Static |
| Kernel setup | `COLLISIONS_HXC.F` | Various | Equilibrium |

### Missing Components

| Feature | Module | Status |
|---------|--------|--------|
| Dielectric computation | `src/pol_function/` | Not called in RT |
| Dynamic screening | - | Missing |
| BSE kernel | `src/bse/` | Not called in RT |
| BSE dynamics | - | Missing |

---

## 10. Physical Implications {#implications}

### What CAN Be Described

1. **Carrier dynamics**: Full occupation evolution
2. **Band-filling bleaching**: State occupation depletion
3. **Hartree screening**: Coulomb repulsion with carriers
4. **Nonlinear current**: Single-particle response to field
5. **Weak-field response**: Perturbative regime

### What CANNOT Be Described

1. **Excitonic bleaching**: Many-body renormalization missing
2. **Screening buildup**: Dynamic polarization not tracked
3. **Plasma effects**: Collective excitations not included
4. **Many-body dephasing**: Beyond mean-field scattering
5. **Strong-field effects**: Breakdown of static kernel approximation

### Quantitative Estimates

When is static kernel valid?

- **Weak field**: $E < 10^5$ V/cm (safely valid)
- **Moderate field**: $10^5$ - $10^7$ V/cm (marginal)
- **Strong field**: $> 10^7$ V/cm (invalid)

---

## 11. Recommended Usage {#recommended_usage}

### Use YAMBO RT For

- Carrier thermalization with realistic scattering
- Transient absorption with band filling
- High-harmonic generation from single-particle dynamics
- Weak- to moderate-field nonlinear optics
- Ultrafast carrier dynamics (ps timescale)

### Do NOT Use For

- Exciton-resolved spectroscopy
- Excitonic bleaching effects
- Strong-field physics requiring dynamic screening
- Many-body correlations beyond mean-field
- Plasmonics with dynamic screening

### Workarounds

1. **For excitonic effects**: Run separate BSE calculations at different fixed carrier densities
2. **For stronger fields**: Reduce field strength and extrapolate
3. **For screening effects**: Use TDDFT mode for better Hartree + semi-local XC
4. **For validation**: Compare with experiment at multiple field strengths

---

## See Also

- [NEGF Theory](YAMBO_RT_NEGF_THEORY.md): Theoretical foundations
- [Screening Theory](YAMBO_RT_SCREENING_THEORY.md): Dielectric response details
- [Kadanoff-Baym Implementation](YAMBO_RT_KADANOFF_BAYM_IMPLEMENTATION.md): Time integration
- [Hamiltonian Construction](YAMBO_RT_HAMILTONIAN_IMPLEMENTATION.md): H(t) details

## References

1. C. Attaccalite et al., *Phys. Rev. B* **84**, 245110 (2011)
2. D. Sangalli et al., *J. Chem. Phys.* **134**, 034115 (2011)
3. A. Marini et al., *Comp. Phys. Comm.* **180**, 1392 (2009)
4. Source code: `src/real_time_drivers/RT_driver.F`, `src/real_time_hamiltonian/RT_Hamiltonian.F`
