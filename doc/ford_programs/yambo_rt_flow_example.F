!
! License-Identifier: GPL
!
! Copyright (C) 2024 The Yambo Team
!
! This is a conceptual program demonstrating the yambo_rt real-time flow
!
program yambo_rt_flow
 !! **Real-Time NEGF Propagation Flow - Conceptual Program**
 !!
 !! This program demonstrates the logical structure and flow of the yambo_rt
 !! real-time Kadanoff-Baym propagation. It shows how the main components
 !! are orchestrated to evolve the non-equilibrium Green's function in time.
 !!
 !! **Main Algorithm**:
 !!   1. **Initialization**: Setup initial Green's function and equilibrium state
 !!   2. **Time Loop**: Iterate over timesteps
 !!      - Compute Hamiltonian H(t) from charge density
 !!      - Extract occupations f_n(k,t)
 !!      - Propagate G(t) → G(t+dt)
 !!      - Apply scattering/collision terms
 !!      - Compute observables (current, polarization, energy)
 !!   3. **Output**: Save trajectories and observables
 !!
 !! See the actual implementation in:
 !!   - `src/real_time_drivers/RT_driver.F` - Main driver
 !!   - `src/real_time_hamiltonian/RT_Hamiltonian.F` - Hamiltonian computation
 !!   - `src/real_time_propagation/RT_Integrator.F` - Time propagation
 !!   - `src/real_time_propagation/RT_occupations_eval.F` - Occupation extraction
 !!   - `src/real_time_control/RT_Observables.F` - Observable computation
 !!
 implicit none
 
 ! Derived types for system description
 ! type(levels)       :: E        !! Electronic structure
 ! type(bz_samp)      :: k        !! K-point grid
 ! type(X_t)          :: X        !! Response functions
 ! type(DIPOLE_t)     :: Dip      !! Dipole elements
 !
 ! ! Time-dependent quantities
 ! complex(SP)        :: G_lesser(:,:,:,:,:)   !! Lesser Green's function G^<(t)
 ! complex(SP)        :: dG_lesser(:,:,:,:,:)  !! G^< derivative
 ! real(SP)           :: f_occ(:,:,:)          !! Occupations f_n(k,t)
 ! complex(SP)        :: H_rt(:,:,:)           !! RT Hamiltonian
 ! type(gauge_field)  :: A_field               !! External field
 
 ! Executable statements
 
 ! STEP 1: Read input and initialize ground state
 ! ================================================
 ! call io_DB1(E, k)           !! Load ground-state database
 ! call setup(E, k)            !! Setup electronic structure
 ! call bz_samp_indexes(E, k)  !! Setup k-point grids
 
 ! STEP 2: Initialize real-time arrays and Green's functions
 ! ===========================================================
 ! call RT_alloc(E, "all")     !! Allocate RT arrays
 ! call RT_initialize(E, k)    !! Initialize G_lesser from Fermi-Dirac
 
 ! STEP 3: Setup external field (laser, THz, etc.)
 ! ================================================
 ! call field_setup(A_field)   !! Configure gauge field
 
 ! STEP 4: Main time integration loop
 ! ===================================
 ! do while (NE_i_time < NE_steps)
 !
 !   ! Update Hamiltonian at time t
 !   call RT_Hamiltonian(dG_lesser, A_field, E, k)
 !   !    Computes: H^RT(t) = H^EQ + ΔΣ^Hxc[G(t)] + U^ext(t)
 !   !    where ΔΣ^Hxc includes dynamic Hartree, XC, and collision corrections
 !
 !   ! Extract occupations from G^<(t,t)
 !   call RT_occupations_eval(dG_lesser)
 !   !    Computes: f_n(k,t) = Im[G^<_nn(k,t,t)]
 !   !    Handles both electrons and holes
 !
 !   ! Propagate: G^<(t) → G^<(t+dt)
 !   call RT_Integrator(G_lesser, dG_lesser, A_field, E, k)
 !   !    Uses predictor-corrector or RK4 scheme
 !   !    Solves Kadanoff-Baym: i∂_t G^< = [H, G^<] + S_xc
 !
 !   ! Apply scattering/collisions (if enabled)
 !   if (l_use_collisions) then
 !      call RT_relaxation(E, k, NE_time)
 !      !   Applies collision kernel: S_xc = K_eq * dG
 !   endif
 !
 !   ! Compute observables at time t
 !   call RT_Observables(E, k, dG_lesser, observables)
 !   !    Computes:
 !   !      - Current density: j(t) = Tr[v·ρ(t)]
 !   !      - Polarization: P(t) = Tr[r·ρ(t)]
 !   !      - Total energy: E_tot(t)
 !   !      - Number of carriers: N(t)
 !
 !   ! Write output and save restart files
 !   call RT_output_and_IO_driver('write', E, k, observables)
 !
 !   ! Update time
 !   NE_i_time = NE_i_time + 1
 !   NE_time = NE_time + RT_dyn_step
 !
 ! enddo
 
 ! STEP 5: Finalize and output
 ! ============================
 ! call RT_output_and_IO_driver('close', E, k, observables)
 ! call RT_free()  !! Deallocate RT arrays
 
 ! STEP 6: Post-processing with ypp (in separate program)
 ! ========================================================
 ! Analyze RT trajectories:
 !   - Compute spectra from time-dependent dipole
 !   - Extract band populations evolution
 !   - Compute high-harmonic generation spectra
 !   - Generate time-resolved plots
 
 ! Program flow summary (text output)
 print *, "=========================================="
 print *, "  Real-Time NEGF Propagation"
 print *, "=========================================="
 print *, ""
 print *, "Main Entry Points:"
 print *, "  yambo_rt:    Real-time propagation driver"
 print *, "  ypp:         Post-processing and analysis"
 print *, ""
 print *, "Key Subroutines:"
 print *, "  RT_driver:              Main RT loop orchestrator"
 print *, "  RT_Hamiltonian:         Compute H(t) from density"
 print *, "  RT_Integrator:          Advance G^<(t) via KBE"
 print *, "  RT_occupations_eval:    Extract f_n(k,t) from G^<"
 print *, "  RT_Observables:         Compute j(t), P(t), E(t)"
 print *, "  RT_relaxation:          Apply collision kernel"
 print *, ""
 print *, "Theory References:"
 print *, "  - Kadanoff-Baym Equations (KBE)"
 print *, "  - Non-Equilibrium Green's Functions (NEGF)"
 print *, "  - Time-Dependent Hartree or DFT"
 print *, ""
 print *, "See doc/markdown/ for complete theory and"
 print *, "implementation documentation."
 print *, "=========================================="
 
end program yambo_rt_flow
