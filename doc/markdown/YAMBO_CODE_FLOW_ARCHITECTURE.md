---
title: YAMBO Code Architecture and Flow Documentation
author: YAMBO Team
date: 2024-12-04
---

# YAMBO Code Architecture & Flow Documentation

## Overview

This document describes the overall flow and architecture of YAMBO, covering the three main entry points: **yambo**, **yambo_rt**, and **ypp**.

---

## Main Entry Points

### 1. YAMBO - Ground State & Response Functions

**File**: `driver/yambo.F`

**Purpose**: Calculate ground state properties and response functions (GW, BSE, optics, etc.)

**Main Flow**:
```
yambo()
├─ SET_defaults()           ← Initialize all configurations
├─ INIT()                   ← Read ground state DB and setup
├─ io_DB1()                 ← Load database
├─ [Calculation Modules]    ← Based on user input:
│  ├─ X_dielectric_matrix() ← Polarization/Screening (GW)
│  ├─ BSE()                 ← Bethe-Salpeter Equation
│  ├─ NL_Optics()           ← Nonlinear optical response
│  ├─ DIPOLE_build()        ← Dipole elements
│  └─ Collisions()          ← Scattering processes
└─ Output Files             ← Save results
```

**Key Decisions**:
- Line 17-20: Check enabled modules (`l_real_time`, `l_gw0`, `l_bse`, etc.)
- Line 83-89: Initialize system with defaults and configuration transfer
- Main conditional branches based on which calculations to run

**Output**: Database files (R_*, o.*, *, %, &) for response functions

---

### 2. YAMBO_RT - Real-Time Propagation

**File**: `src/real_time_drivers/RT_driver.F`

**Purpose**: Propagate electron density matrix in real time under external fields

**Main Flow**:
```
RT_driver()
├─ RT_initialize()          ← Setup initial Green's function and Hamiltonian
│  ├─ RT_G_lesser_init()    ← Initialize G^< from Fermi-Dirac
│  ├─ RT_occupations_update() ← Extract occupations
│  └─ RT_H_setup()          ← Build initial Hamiltonian
│
├─ RT_time_integration_loop()
│  ├─ RT_Hamiltonian()      ← Compute H(t) at each step
│  │  ├─ el_density_matrix() ← Extract ρ(r,t) from G(t)
│  │  ├─ V_Hartree()        ← Solve Poisson: V^H = ∫ ρ(r')/|r-r'|
│  │  ├─ XC_potential()     ← Compute V^xc (TDDFT or static)
│  │  └─ RT_apply_field()   ← Add external perturbation
│  │
│  ├─ Propagator()          ← Advance G(t) → G(t+dt)
│  │  ├─ PREDIC_CORREC()    ← Predictor-corrector step
│  │  └─ Collisions_RT()    ← Apply scattering kernel
│  │
│  ├─ RT_occupations_eval() ← Extract occupations at t+dt
│  │  └─ f_n(k,t) = Im[G^<_nn(k,t,t)]
│  │
│  └─ RT_observables()      ← Compute output quantities
│     ├─ Current           ← j(t)
│     ├─ Polarization      ← P(t)
│     └─ Energy            ← E(t)
│
└─ Output: Time-dependent arrays
   ├─ Occupations vs time
   ├─ Current vs time
   ├─ Polarization vs time
   └─ Restart files
```

**Key Decisions**:
- Hartree potential: updated dynamically from charge density
- XC potential: fully dynamic (TDDFT) or static (collision kernel mode)
- Screening: Uses pre-computed GW kernel (NOT updated in RT)
- Collision effects: Applied as static kernel × ΔG

**Coupling to Hamiltonian**:
$$i\hbar\frac{\partial G^<}{\partial t} = [H(t), G^<(t)] + \text{collisions}$$

**Output**: RT trajectory files with observables

---

### 3. YPP - Post-Processing & Analysis

**File**: `driver/ypp.F`

**Purpose**: Analyze and post-process results from yambo and yambo_rt calculations

**Main Flow**:
```
ypp()
├─ SET_defaults()           ← Initialize YPP configurations
├─ Load databases           ← io_DB1(), io_QINDX(), io_GROT()
├─ [Analysis Modules]       ← Based on user flags:
│  ├─ YPP_RealTime()        ← Analyze RT trajectories (if compiled with _YPP_RT)
│  │  ├─ Read RT output files
│  │  ├─ Extract observables (current, polarization, etc.)
│  │  ├─ Fourier analysis → Spectrum
│  │  └─ Phase analysis
│  │
│  ├─ YPP_Electrons()       ← Band structure analysis
│  │  ├─ Band diagram plotting
│  │  └─ DOS calculation
│  │
│  ├─ YPP_Excitons()        ← Excitonic spectroscopy
│  │  ├─ Read BSE eigenvectors
│  │  ├─ Build absorption spectrum
│  │  └─ Excitonic wave functions
│  │
│  ├─ YPP_Wannier()         ← Wannier function analysis
│  └─ [Other modules]       ← QP DBs, ELPH, dipoles, etc.
│
└─ Output: Analysis files
   ├─ Spectra (absorption, etc.)
   ├─ Formatted datasets
   └─ Graphical data
```

**Key Decision Points**:
- Line 28-34: Check which YPP modules are compiled
- YPP_RealTime branch: Only active if `_YPP_RT` is defined
- Database manipulation: Create/expand/modify QP databases

**Workflow Integration**: YPP can run standalone or chains with yambo results

---

## Critical Subroutines

### RT_Hamiltonian - Heart of Real-Time
**File**: `src/real_time_hamiltonian/RT_Hamiltonian.F`

**Purpose**: Compute time-dependent Hamiltonian for KBE

**Formula**:
$$H^{RT}(t) = H^{EQ} + \Delta V^{H}(t) + \Delta V^{xc}(t) + U^{ext}(t)$$

**Components**:
1. **Hartree**: $\Delta V^{H}(t) = \int \frac{\Delta\rho(r',t)}{|r-r'|} d^3r'$
2. **XC Potential**: 
   - TDDFT mode: $V^{xc}[n(r,t)]$ (dynamic)
   - Collision kernel: Static $K_{eq}$ applied to $\Delta G(t)$
3. **External Field**: $U^{ext}(t) = E(t) \cdot r$ or velocity gauge

**Flow**:
```
RT_Hamiltonian(dG, A_field) → Ho_plus_Sigma
├─ Extract charge density: ρ(r,t) from dG
├─ Solve Poisson equation: V^H(ρ)
├─ Compute XC potential: V^xc(ρ) if TDDFT
├─ Apply collision kernel: K_eq * dG if SEX/COHSEX
└─ Add external field: A(t) coupling
```

**Key Limitation**: Screening NOT updated (uses equilibrium kernel)

---

### RT_occupations_eval - Occupation Dynamics
**File**: `src/real_time_propagation/RT_occupations_eval.F`

**Purpose**: Extract time-dependent band occupations

**Formula**:
$$f_n(k,t) = \text{Im}[G^<_{nn}(k,t,t)]$$

**Captures**:
- Band-filling (depletion/filling of states)
- Carrier generation from excitation
- Thermal equilibration (if scattering enabled)

---

### Collisions_RT - Scattering Effects
**File**: `src/collisions/COLLISIONS_compose_rt.F`

**Purpose**: Apply collision kernel for relaxation and dephasing

**Self-Energy**:
$$\Sigma^{xc}(t) = K_{eq}[G(t) - G_{ref}]$$

**Where**:
- $K_{eq}$ = equilibrium kernel from GW/BSE pre-computation
- $G(t) - G_{ref}$ = deviation from reference (ground state)
- Result: Self-consistent mean-field scattering

**Status**: STATIC kernel (does not evolve with carrier dynamics)

---

## Data Flow Architecture

### Ground State → Real-Time

```
                    ┌──────────────────┐
                    │   Ground State   │
                    │   Calculation    │
                    │    (yambo)       │
                    └────────┬─────────┘
                             │
                             ↓
                    ┌──────────────────┐
                    │   Save DB1       │
                    │  (k-points,      │
                    │   energy bands)  │
                    └────────┬─────────┘
                             │
                    ┌────────v──────────┐
                    │   GW Screening    │
                    │   BSE Kernels     │
                    │  (if needed)      │
                    └────────┬─────────┘
                             │
                    ┌────────v──────────┐
                    │  Real-Time (RT)   │
                    │  Propagation      │
                    │  (yambo_rt)       │
                    └────────┬─────────┘
                             │
                    ┌────────v──────────┐
                    │   RT Observables  │
                    │  (j, P, E, etc.)  │
                    └────────┬─────────┘
                             │
                    ┌────────v──────────┐
                    │   Analysis        │
                    │   (ypp)           │
                    │   Spectroscopy    │
                    └──────────────────┘
```

### Real-Time Loop (Per Timestep)

```
        G(t) ──────────────┐
         │                 │
         ├─→ Extract ρ ────┤
         │                 │
         ├─→ Hartree ──────┤→ H(t)
         │                 │
         ├─→ XC ───────────┤
         │                 │
         └─→ Field ────────┘
                 │
                 ↓
          Propagator
          ├─ Predictor
          └─ Corrector
                 │
                 ↓
             G(t+dt)
                 │
                 ├→ Occupations
                 ├→ Current
                 ├→ Polarization
                 └→ Observables
```

---

## Module Dependencies

### Essential RT Modules

| Module | Purpose | Dependency |
|--------|---------|-----------|
| `real_time` | RT control, arrays, time integration | Core |
| `real_time_hamiltonian` | H(t) computation | real_time |
| `real_time_propagation` | Green's function propagation | real_time |
| `collision_ext` | Collision kernel setup | GW pre-computed |
| `fields` | External field control | real_time |
| `FFT_m` | Real-space Poisson solver | Core |
| `electrons` | Occupation management | Core |

### Input/Output

| Module | Purpose |
|--------|---------|
| `IO_m` | Database I/O |
| `IO_int` | Internal control |
| `RT_control` | RT I/O settings |

---

## Key Physics Features

### Implemented ✓
- [x] Occupation evolution (band-filling, excitation)
- [x] Hartree potential (dynamically updated)
- [x] External field coupling (velocity/length gauge)
- [x] Collisions (static kernel from GW)
- [x] Time integration (predictor-corrector)

### Limitations ⚠
- [ ] Dynamic screening (uses static $\varepsilon$)
- [ ] Excitonic effects (no dynamic BSE kernel)
- [ ] Band-gap renormalization (static approximation)

### Physical Domain
- **Valid for**: Weak to moderate pump fields
- **Screening approximation**: Equilibrium or time-dependent Hartree-Fock
- **Scattering**: Simplified collision kernels
- **Time scale**: Femtoseconds to picoseconds

---

## Execution Paths

### Path 1: Static Response (GW, BSE, Optics)
```
yambo
├─ INIT (ground state setup)
├─ X_dielectric_matrix (RPA screening)
├─ Correlation effects (if enabled)
└─ Save response databases
```

### Path 2: Real-Time Dynamics
```
yambo_rt
├─ Read ground state + GW kernel
├─ RT_initialize (prepare G and H)
├─ for each timestep:
│  ├─ RT_Hamiltonian (build H(t))
│  ├─ Propagator (advance G)
│  └─ Observables (extract data)
└─ Save RT trajectory
```

### Path 3: Post-Processing
```
ypp
├─ Load yambo/yambo_rt outputs
├─ Analyze/transform data
└─ Export formatted results
```

---

## File Organization

```
Lumen/
├── driver/
│   ├── yambo.F              ← Main ground-state driver
│   └── ypp.F                ← Post-processing driver
│
├── src/
│   ├── real_time_drivers/
│   │   └── RT_driver.F      ← RT loop orchestration
│   │
│   ├── real_time_hamiltonian/
│   │   └── RT_Hamiltonian.F ← H(t) computation
│   │
│   ├── real_time_propagation/
│   │   ├── RT_occupations_eval.F  ← f_n(k,t)
│   │   └── Propagator.F           ← G(t) advance
│   │
│   ├── collisions/
│   │   └── COLLISIONS_compose_rt.F ← Scattering
│   │
│   ├── response_function/
│   │   └── X_dielectric_matrix.F   ← Screening (GW)
│   │
│   ├── bse/
│   │   └── K_*.F             ← BSE kernels
│   │
│   └── ... (other modules)
│
├── ypp/
│   └── src/
│       ├── YPP_RealTime*.F   ← RT analysis
│       ├── YPP_*.F           ← Other YPP modules
│       └── ...
│
└── include/
    └── y_*.h                ← Header files
```

---

## Control Flow: A Typical Simulation

1. **Setup** (yambo or yambo_rt executable)
   - Parse input file
   - Initialize MPI parallelization
   - Call main driver function

2. **Ground State** (if needed, yambo)
   - Load DFT calculation
   - Compute response functions
   - Save databases

3. **Real-Time** (yambo_rt)
   - Load ground state + screening
   - Initialize G^< and Hamiltonian
   - Timestep loop:
     - Compute H(t)
     - Propagate G(t)
     - Extract observables
   - Save trajectory

4. **Analysis** (ypp)
   - Load RT output
   - Fourier transform (j → conductivity, σ)
   - Build spectra (absorption, etc.)
   - Export results

---

## Documentation References

- **Theory**: [YAMBO_RT_NEGF_THEORY.md](YAMBO_RT_NEGF_THEORY.md)
- **Implementation**: [YAMBO_RT_NEGF_IMPLEMENTATION.md](YAMBO_RT_NEGF_IMPLEMENTATION.md)
- **Screening**: [YAMBO_RT_SCREENING_THEORY.md](YAMBO_RT_SCREENING_THEORY.md)
- **Quick Reference**: [YAMBO_RT_NEGF_QUICK_REFERENCE.md](YAMBO_RT_NEGF_QUICK_REFERENCE.md)

---

**Last Updated**: 2024-12-04  
**Status**: Complete - Ready for FORD integration
