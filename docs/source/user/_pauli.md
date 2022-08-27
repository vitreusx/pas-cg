# Pauli exclusion FF

The Pauli exclusion force field is described by:

$$ \begin{align*}V &= \sum \mathsf{Pauli}(r)\\ \mathsf{Pauli}(r)
&=\begin{cases}\varepsilon\left[\left(\frac{r_{\mathrm{excl}}}{r}\right)^{12}-2\left(\frac{r_{\mathrm{excl}}}{r}\right)^{6}+1\right]
& r \leq r_{\mathrm{excl}}\\ 0 & \mathrm{otherwise} \end{cases}\end{align*} $$

The summation is goes over the non-local pairs which are not in contact - if enabled, the QA and PID potentials handle too-close-residues in their own fashion.

## Parameter file entry

```yaml
# Parameters for the Pauli exclusion force.
Pauli exclusion:
  # Whether it's enabled.
  enabled: true

  # Depth of the force.
  depth: 1.0 eps
```