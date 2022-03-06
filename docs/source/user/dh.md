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
constant Debye-Hueckel:
  enabled: boolean
  screening distance: quantity [L]
  permittivity: quantity [(dim of \varepsilon_0)]

relative Debye-Hueckel:
  enabled: boolean
  screening distance: quantity [L]
  factor: quantity [(dim of \varepsilon_0)/L]
```

Only one of these can be enabled at the same time.