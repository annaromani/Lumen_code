---
title: Documentation Pages
---

# Documentation Pages {#doc_index}

Documentation organized by **Theory → Code → Visual Guides** for logical learning progression.

---

## THEORY & FUNDAMENTALS

Theoretical foundations and mathematical frameworks to understand LUMEN's physics.

### Equilibrium Green's Functions

- [YAMBO Theory Documentation](YAMBO_THEORY_DOCUMENTATION.html) — Foundational concepts and mathematical formulations
- [Theory Code Examples](YAMBO_THEORY_CODE_EXAMPLES.html) — Practical code implementations of theory concepts

### Real-Time & Non-Equilibrium Systems

- [Real-Time TDDFT Theory](YAMBO_RT_THEORY_DOCUMENTATION.html) — Time-dependent density functional theory fundamentals
- [NEGF Theory Guide](YAMBO_RT_NEGF_THEORY.html) — Non-equilibrium Green's functions fundamentals
- [Screening in RT Systems](YAMBO_RT_SCREENING_THEORY.html) — Screening effects in real-time simulations
- [Kadanoff-Baym Theory & Implementation](KADANOFF_BAYM_IMPLEMENTATION.html) — KBE derivation and code mapping

### Exciton-Phonon Coupling

- [Exciton-Phonon Theory](EXCPH_THEORY.html) — Complete theoretical framework for exciton-phonon interactions
- [Phonon-Assisted Transitions](EXCPH_IMPLEMENTATION.html) — Theory of phonon-assisted optical processes
- [Exciton Lifetimes](EXCPH_IMPLEMENTATION.html) — Non-radiative decay rates from phonon scattering

### Quick References

- [Formula Quick Reference](YAMBO_FORMULAS_QUICK_REFERENCE.html) — Essential equations at a glance
- [NEGF Implementation Quick Reference](YAMBO_RT_NEGF_QUICK_REFERENCE.html) — Quick summary of NEGF features

---

## IMPLEMENTATION & ARCHITECTURE

Code structure, algorithms, and implementation details for developers.

### Exciton-Phonon Coupling

- [Exciton-Phonon Implementation](EXCPH_IMPLEMENTATION.html) — Computational methods and algorithms for coupling calculation
- [Matrix Element Computation](EXCPH_IMPLEMENTATION.html) — How gkkp matrix elements are calculated
- [Optical Properties Calculation](EXCPH_IMPLEMENTATION.html) — Computing absorption and luminescence spectra
- [Lifetime Calculations](EXCPH_IMPLEMENTATION.html) — Non-radiative decay rates from scattering

### Real-Time Simulations

- [RT NEGF Implementation Analysis](YAMBO_RT_NEGF_IMPLEMENTATION.html) — What features are actually implemented
- [RT Screening Implementation](YAMBO_RT_SCREENING_IMPLEMENTATION.html) — Screening algorithms in code
- [Two-Time Green's Function Reconstruction](TWO_TIME_GREEN_FUNCTION_RECONSTRUCTION.html) — Advanced reconstruction techniques
- [Two-Time Summary](TWO_TIME_RECONSTRUCTION_SUMMARY.html) — Quick summary of two-time approach

### System Architecture

- [Code Flow & Architecture](YAMBO_CODE_FLOW_ARCHITECTURE.html) — Main entry points and control flow
- [I/O System Architecture](IO_System_Architecture.html) — File I/O and database organization
- [I/O Quick Reference](IO_Quick_Reference.html) — I/O API and file formats

### Memory & Allocation

- [Memory Allocation System](YAMBO_ALLOC_DOCUMENTATION.html) — Dynamic memory management design
- [Allocation Architecture](YAMBO_ALLOC_ARCHITECTURE.html) — Allocator internals and patterns
- [Allocation Examples](YAMBO_ALLOC_EXAMPLES.html) — Code examples for memory allocation
- [Allocation Quick Reference](YAMBO_ALLOC_QUICK_REFERENCE.html) — API quick guide

### User Configuration

- [Input Parameters Reference](INPUT_FLAGS.html) — All ~177 input flags documented

---

## VISUAL GUIDES & FLOWCHARTS

Diagrams, flowcharts, and visual learning materials.

### Real-Time Dynamics

- [Real-Time Propagation Visual Guide](YAMBO_RT_VISUAL_GUIDE.html) — Step-by-step RT algorithm visualization
- [YPP Real-Time Flowcharts](YPP_RT_FLOWCHARTS.html) — Post-processing workflows and pipelines
- [YPP Real-Time Documentation](YPP_RT_DOCUMENTATION.html) — Complete YPP RT module guide

### Theory Visualization

- [Theory Visual Guide](YAMBO_THEORY_VISUAL_GUIDE.html) — Diagrams and illustrations of key concepts

### System Diagrams

- [I/O System Flowchart](IO_System_Flowchart.html) — Visual I/O control flow

---

## Main Documentation

For comprehensive overview and getting started guide, see [Yambo Documentation](mainpage.html)
