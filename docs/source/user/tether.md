# Tether pair FF

The tether pair force field is described by:

$$V = \sum
\frac{1}{2}\left[\mathtt{H1}(r_{12} - d)^2+\mathtt{H2}(r_{12}-d)^4\right]$$

where the summation is over the tethers of the model and the value of $d$ is
either the native value $d_0$, if such is present, or $d_{\textrm{def}}$
otherwise.

## Parameter file entry

```yaml
tether forces:
  enabled: boolean
  H1: quantity [E/L^2]
  H2: quantity [E/L^2]
  default bond length: quantity [L]
```