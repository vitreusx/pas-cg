# Pauli exclusion FF

The Pauli exclusion force field is described by:

$$ \begin{align*}V &= \sum \mathsf{Pauli}(r)\\ \mathsf{Pauli}(r)
&=\begin{cases}\varepsilon\left[\left(\frac{r_{\mathrm{excl}}}{r}\right)^{12}-2\left(\frac{r_{\mathrm{excl}}}{r}\right)^{6}+1\right]
& r \leq r_{\mathrm{excl}}\\ 0 & \mathrm{otherwise} \end{cases}\end{align*} $$

The summation is over **all** non-local residue pairs, irrespective of whether
they are in (QA) contact or not.

## Parameter file entry

```yaml
Pauli exclusion:
  enabled: boolean
  exclusion radius: quantity [L]
  depth: quantity [E]
```