# Reaction Term and Its YAML Definition

This page explains the reaction term $R_k$ used in the model equations and then maps that formulation to the YAML reactions block.

## 1. Reaction-Side Definition First

In the component mass balance, $R_k(\mathbf{C}, x, t)$ is the local net reaction source or sink for species $k$.

For a reaction network with reactions $r=1,\dots,N_R$:

$$
R_k = \sum_{r=1}^{N_R} \nu_{k,r}\cdot \rho_r
$$

where:

- $\rho_r$ is the net rate of reaction $r$
- $\nu_{k,r}$ is the stoichiometric coefficient of species $k$ in reaction $r$

Sign convention:

- $\nu_{k,r} < 0$: species $k$ is consumed in reaction $r$
- $\nu_{k,r} > 0$: species $k$ is produced in reaction $r$

## 2. Worked Equation Example (from equilibrium.yaml)

Use reaction 13 from equilibrium.yaml:

$$
\mathrm{H_2S} \rightleftharpoons \mathrm{HS^-} + \mathrm{H^+}
$$

One possible forward-backward representation is:

For a reversible elementary reaction, this comes from the law of mass action. In general, for
$A \rightleftharpoons B + C$, the net rate is written as

$$
\rho = \rho_{+} - \rho_{-} = k_{+}\cdot a_A - k_{-}\cdot a_B\cdot a_C
$$

with activities (or concentrations in the simplified model form) $a_i$ and two kinetic constants.
The idea is: one term drives the reaction in one direction, the other term drives it back, and
their difference is the observable net rate.

Applied to the equilibrium.yaml entry, the helper terms `F13` and `B13` are exactly these two
directional contributions. Their names are just labels; what matters is the consistent use in
`rate: "F13 - B13"` together with stoichiometry and the equilibrium constraint.

$$
\rho_{13} = F13 - B13 = kfs\cdot hplus\cdot hs - kbs\cdot h2s
$$

Then the reaction contribution to each species term is:

$$
R_{h2s}^{(13)} = +1\cdot\rho_{13}, \qquad
R_{hs}^{(13)} = -1\cdot\rho_{13}, \qquad
R_{hplus}^{(13)} = -1\cdot\rho_{13}
$$

For equilibrium handling, the same reaction can be written as an algebraic constraint:

$$
hplus\cdot hs - keqs\cdot h2s = 0
$$

This is directly linked to the forward-backward form: at equilibrium, net rate vanishes
($\rho_{13}=0$), so $F13=B13$. Rearranging gives the same concentration relationship encoded by
`equilibrium_constraint`.

## 3. Mapping the Equation to YAML

The mathematical objects map directly to YAML fields:

- $\rho_r$ maps to rate
- helper terms (such as $F13$, $B13$) map to rate_components
- $\nu_{k,r}$ maps to stoichiometry coefficients
- equilibrium algebraic condition maps to equilibrium_constraint

Example entry:

```yaml
  - id: 13
    name: "sulfide_dissociation"
    description: "H2S <-> HS- + H+"
    equilibrium: true
    rate: "F13 - B13"
    rate_components:
      F13: "kfs * hplus * hs"
      B13: "kbs * h2s"
    equilibrium_constraint: "hplus * hs - keqs * h2s"
    stoichiometry:
      h2s: 1
      hs: -1
      hplus: -1
```

## 4. Typical Rate-Law Patterns

| Type | Pattern for rate | Example |
| --- | --- | --- |
| Elementary mass action | `k * A * B` | `k14*o2*(hs + h2s)` |
| Saturation / Monod-like | `k * S/(K + S)` | `kfox * ch2o * (o2/kmo2)` |
| Inhibited rate | `k * f_substrate * (1 - f_inhib)` | `hso4*(1-fo2-fno3-fmn4-ffe3)*(...)` |
| Forward-backward form | `F - B` | `F13 - B13` |
| Precipitation/dissolution | `k * (omega - 1)` | `(k17_1*sw17 + k17_2*mnco3*(1-sw17))*(omega_mn-1)` |

## 5. Minimal Checklist for YAML reactions

- Each reaction has id, name, rate, and stoichiometry.
- Stoichiometric signs are consistent with the chosen rate direction.
- Every symbol in rate is defined as a species, parameter, or rate component.
- If equilibrium is true, equilibrium_constraint is present.
- Units and scaling terms (for example SD) are consistent across the model.