---
title: Screening and Dielectric Response - Theory
---

# Screening and Dielectric Response - Theory {#rt_screening_theory}

**Objective**: Understand screening effects and dielectric response in real-time simulations.

**Prerequisites**:
- [NEGF Theory](YAMBO_RT_NEGF_THEORY.md)
- [Kadanoff-Baym Equations](YAMBO_RT_KADANOFF_BAYM_THEORY.md)

**Reading Time**: 25-35 minutes

---

## 1. Overview {#overview}

Screening is the reduction of the Coulomb interaction between charged particles due to the polarization response of the medium. In real-time dynamics, screening adapts as the carrier population changes, affecting:

- Hartree repulsion strength
- Exchange-correlation kernel
- Exciton binding energies  
- Plasma oscillations

Understanding the difference between static and dynamic screening is crucial for real-time simulations.

---

## 2. The Coulomb Interaction and Screening {#coulomb_screening}

### Bare Coulomb Interaction

The bare Coulomb kernel is:

$$v(q) = \frac{4\pi e^2}{q^2}$$

where $q$ is the wavevector. This diverges as $q \to 0$ and represents unscreened Coulomb interaction.

### Response Function and Polarization

When charges respond to an external field, they create a polarization:

$$\mathbf{P}(\mathbf{r},t) = -\int d^3r' \chi(\mathbf{r},\mathbf{r}',t) E(\mathbf{r}',t)$$

where $\chi$ is the polarization response tensor.

In reciprocal space and frequency domain:

$$P(q,\omega) = -\chi(q,\omega) E(q,\omega)$$

### Lindhard Polarization

The fundamental response comes from particle-hole excitations. The bare polarization function is:

$$\chi_0(q,\omega) = -\sum_{nk} \frac{f_n(k) - f_m(k+q)}{E_m(k+q) - E_n(k) - \omega}$$

where the sum is over occupied ($n$) and empty ($m$) states. This formula applies in equilibrium where occupations are time-independent.

---

## 3. Dielectric Function and Screening {#dielectric_function}

### Definition of Dielectric Function

The dielectric function relates the total electric field to the applied field:

$$\mathbf{E}(\mathbf{r},\omega) = \varepsilon(\omega) E_0(\mathbf{r},\omega)$$

In reciprocal space:

$$\varepsilon(q,\omega) = 1 - v(q)\chi(q,\omega)$$

where $\chi$ includes all correlation effects.

### Screened Coulomb Interaction

The effective Coulomb interaction between particles is:

$$W(q,\omega) = \frac{v(q)}{\varepsilon(q,\omega)}$$

This is the **screened** interaction that appears in the self-energy and BSE kernel.

### Static vs. Dynamic Screening

**Static screening**:
$$\varepsilon_0(q) = \varepsilon(q,\omega=0) = 1 - v(q)\chi_0(q,0)$$

Used when interactions are instantaneous or when frequency dependence is negligible.

**Dynamic screening**:
$$\varepsilon(q,\omega) = 1 - v(q)\chi(q,\omega)$$

Full frequency-dependent response, important for spectral properties.

---

## 4. Real-Time Screening: Time Dependence {#real_time_screening}

### Time-Dependent Lindhard Function

In real-time, occupations evolve: $f_n(k,t)$. The polarization becomes:

$$\chi(q,\omega,t) = -\sum_{nk} \frac{f_n(k,t) - f_m(k+q,t)}{E_m(k+q) - E_n(k) - \omega}$$

The response function now depends explicitly on time through the occupations.

### Dynamic Dielectric Function

$$\varepsilon(q,\omega,t) = 1 - v(q)\chi(q,\omega,t)$$

Screening adapts to the changing carrier distribution.

### Physical Meaning

As carriers are photocreated:

1. **Increase in polarization**: More electrons and holes to screen external fields
2. **Plasma frequency shift**: $\omega_p^2(t) \propto n_{\text{carriers}}(t)$
3. **Screening buildup**: Occurs on the timescale of carrier dynamics
4. **Band gap renormalization**: Due to many-body screening

---

## 5. Hartree Screening {#hartree_screening}

### Classical Treatment

The Hartree potential from electron density $\rho(r,t)$ is:

$$V^H(r,t) = \int \frac{\rho(r',t)}{|{\mathbf{r}} - {\mathbf{r}}'|} d^3r'$$

This is exact for classical electrostatics with instantaneous response (infinite dielectric constant).

### In Fourier Space

$$V^H(q,t) = v(q) \rho(q,t) = \frac{4\pi}{q^2} \rho(q,t)$$

The Hartree potential screening comes from the density response:

- Higher carrier density → stronger Hartree repulsion
- Dynamic adaptation as $n(t)$ changes

### Hartree Mean-Field

In mean-field theory:

$$\Delta E^H_n(t) = \int V^H(r,t) |{\psi_n(r)}|^2 d^3r$$

includes density-dependent energy shifts.

---

## 6. Hartree-Fock Exchange Screening {#hf_screening}

### Exchange Self-Energy

The exchange self-energy is:

$$\Sigma^x_{\text{ij}}(1,2) = -\sum_k \langle ik|jk \rangle \frac{1}{|{\mathbf{r}}_1 - {\mathbf{r}}_2|}$$

In equilibrium, this depends on occupied states only.

### Screening in Exchange

Exchange is not directly screened by the dielectric response (no frequency dependence in bare exchange). However, the effective exchange kernel can be written as:

$$K^x_{\text{ij},kl} = \langle ik|v|jl \rangle - \langle il|v|jk \rangle$$

where the Coulomb interaction $v$ is the bare interaction.

### Self-Consistent Screening

The exchange-correlation self-energy can be written more generally as:

$$\Sigma^{\text{xc}} = \Sigma^x + \Sigma^c$$

where correlation self-energy involves screening.

---

## 7. Screened Exchange (SEX) and GW {#sex_gw}

### Screened Coulomb in Self-Energy

The more accurate self-energy uses the **screened** Coulomb interaction:

$$\Sigma^{\text{SEX}}_{\text{ij}} = -\sum_k \langle ik|W|jk \rangle f_k$$

where $W = v/\varepsilon$ is the screened interaction at $\omega=0$ (static approximation).

### GW Approximation

The full GW self-energy includes dynamic screening:

$$\Sigma^{\text{GW}}(i,j,\omega) = \sum_k \langle ik|W(\omega)|jk \rangle$$

where $W(\omega)$ is frequency-dependent:

$$W(q,\omega) = \frac{v(q)}{\varepsilon(q,\omega)}$$

### Static vs. Dynamic in GW

**Static GW** (COHSEX):
- Use $W(q,\omega=0)$ only
- Fast to compute
- Reasonable approximation for excitations

**Dynamic GW**:
- Full frequency dependence
- Expensive but accurate
- Captures spectral features

---

## 8. Vertex Effects {#vertex_effects}

### Beyond RPA

The random phase approximation (RPA) includes only bubble diagrams:

$$\chi_0^{\text{RPA}}(q,\omega) = \sum_{nk} \frac{f_n - f_m}{E_n - E_m - \omega}$$

Higher-order effects come from vertex corrections.

### Vertex Function

The vertex function $\Gamma$ modifies the polarization:

$$\chi(q,\omega) = \chi_0(q,\omega) + \chi_0 \Gamma \chi$$

This includes ladder diagrams (repeated interactions) and other non-RPA effects.

### Physical Consequences

Vertex corrections modify:
- Screening strength
- Excitonic binding energies
- Non-linear optical response

In real-time dynamics, these are typically approximated or neglected.

---

## 9. Plasma Oscillations {#plasma}

### Drude Model

For free carriers with collision frequency $\nu$:

$$\varepsilon(\omega) = 1 - \frac{\omega_p^2}{\omega^2 + i\nu\omega}$$

where the plasma frequency is:

$$\omega_p^2 = \frac{n e^2}{m \varepsilon_0}$$

### Plasma Frequency with Density

For time-dependent carrier density:

$$\omega_p^2(t) = \frac{n_{\text{carriers}}(t) e^2}{m \varepsilon_0}$$

As carriers are photoexcited, the plasma frequency increases, leading to:
- Larger screening at lower frequencies
- Shift of collective modes
- Enhanced polarizability

### Collective Excitations

Collective oscillations of the electron gas occur near $\omega \approx \omega_p$. These are important for:
- Plasmonic effects
- Optical response
- Terahertz dynamics

---

## 10. Excitonic Screening {#excitonic_screening}

### Exciton Definition

An exciton is a bound state of an electron-hole pair, with binding energy:

$$E_{\text{b}} = E_g - E_{\text{ex}}$$

where $E_g$ is the band gap and $E_{\text{ex}}$ is the exciton energy.

### Screening of Exciton Binding Energy

The binding energy depends on screening:

$$E_{\text{b}} \propto \frac{1}{\varepsilon^2}$$

Higher screening (larger $\varepsilon$) reduces binding energy.

### Dynamic Screening of Excitons

As carriers are created:
1. Polarization increases
2. Dielectric constant $\varepsilon$ increases
3. Binding energy decreases (exciton is **ionized**)
4. Exciton is converted to free carriers

This is the **screening-induced exciton dissociation**.

### Band Gap Renormalization

Carrier-induced screening also shifts the band gap:

$$\Delta E_g(n_c) = -\alpha \sqrt{n_c}$$

where $\alpha$ is a material-dependent constant and $n_c$ is carrier density.

---

## 11. Comparison: Static vs. Dynamic Screening {#comparison}

### Static Screening Approximation

Use $\varepsilon(q,\omega=0)$ or $\varepsilon(q)$ without frequency dependence.

**Assumptions**:
- Response is instantaneous
- Frequency-dependent effects negligible
- Valid for low frequencies and long wavelengths

**Advantages**:
- Computationally efficient
- Good for equilibrium properties
- Standard in GW/BSE

**Limitations**:
- Misses frequency dependence
- Invalid for high-frequency processes
- Poor for spectral functions

### Full Dynamic Screening

Include frequency dependence: $\varepsilon(q,\omega,t)$

**Advantages**:
- Captures spectral features
- Valid over broad frequency range
- Accurate for dynamics

**Limitations**:
- Computationally very expensive
- Requires memory in time/frequency
- Complex implementation

### Time-Dependent Screening

Update screening as occupation changes: $\varepsilon(q,\omega,t)$

**Physical picture**:
- Screening adapts on carrier thermalization timescale
- Plasma frequency shifts with density
- Band gap renormalization included

**Computational cost**:
- Recompute $\chi(t)$ at each step
- FFT scaling in frequency domain
- Prohibitive for large systems

---

## 12. Applications and Regimes {#applications}

### When Static Screening is Valid

1. **Weak perturbations**: $\Delta n_e \ll n_e$
2. **Equilibrium properties**: Band structures, bandgaps
3. **Low carrier density**: Isolated excitations
4. **Long-wavelength limit**: $qa \ll 1$

### When Dynamic Screening is Needed

1. **Spectral functions**: Requires frequency resolution
2. **Core excitons**: Strong correlations
3. **Multiple excitations**: Frequency-dependent screening crucial
4. **Non-linear processes**: Harmonic generation

### When Time-Dependent Screening is Important

1. **Ultrafast carrier dynamics**: fs timescale
2. **Plasmonic effects**: Collective mode evolution
3. **Pump-probe spectroscopy**: Time-resolved screening
4. **Strong-field ionization**: Rapid density changes

---

## See Also

- [NEGF Theory](YAMBO_RT_NEGF_THEORY.md): Complete NEGF framework
- [Screening Implementation](YAMBO_RT_SCREENING_IMPLEMENTATION.md): How screening is computed
- [NEGF Implementation](YAMBO_RT_NEGF_IMPLEMENTATION.md): Real-time implementation status

## References

1. H. Haug and A.-P. Jauho, *Quantum Kinetics in Transport and Optics of Semiconductors* (Springer, 2008)
2. L.V. Keldysh, *Sov. Phys. JETP* **20**, 1018 (1965)
3. E.M. Lifshitz and L.P. Pitaevskii, *Physical Kinetics* (Pergamon Press, 1981)
4. G. Mahan, *Many-Particle Physics* (Plenum Press, 1990)
5. C. Attaccalite et al., *Phys. Rev. B* **84**, 245110 (2011)
