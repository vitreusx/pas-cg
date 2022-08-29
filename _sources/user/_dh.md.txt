# Debye-Hueckel FF

The Debye-Hueckel force field exists in two variants: with _constant_ and
*relative* electric permittivity of the medium. Both are described by:

$$ V = \sum \frac{q_1 q_2}{4\pi \varepsilon(r)} \frac{\exp(-r/s)}{r} $$

where the summation is over all non-local residue pairs. The variants differ in
the choice of the $\varepsilon(r)$ function.

## Constant DH force field

For the constant variant, the choice of $\varepsilon(r)$ is:

$$ \varepsilon(r) = \varepsilon_\mathrm{const} $$

## Relative DH force field

For the relative variant, the choice of $\varepsilon(r)$ is:

$$ \varepsilon(r) = A_\varepsilon r $$

## Parameter file entry

```yaml
# Parameters for the Debye-Hueckel (electrostatics) potential.
Debye-Hueckel:
  # Whether it's enabled.
  enabled: true
  # Variant of the D-H potential. Options:
  # - "constant", in which permittivity is constant;
  # - "relative", in which permittivity is described as eps = A r, with "A" a
  #   permittivity factor.
  variant: constant

  # Screening distance for the D-H potential.
  screening distance: 15.0 A

  # A special cutoff distance for the electrostatics.
  cutoff distance: 20.0 A

  # Parameters for the constant-permittivity version.
  constant variant params:
    # Permittivity of the medium.
    permittivity: 80.0 eps_0

  # Parameters for the relative-permittivity version.
  relative variant params:
    # Permittivity factor, as described above.
    permittivity factor: 4.0 eps_0/A
```