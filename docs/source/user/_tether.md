# Tether pair FF

The tether pair force field is described by:

$$V = \sum
\frac{1}{2}\left[\mathtt{H1}(r_{12} - d)^2+\mathtt{H2}(r_{12}-d)^4\right]$$

where the summation is over the tethers of the model and the value of $d$ is
either the native value $d_0$, if such is present, or $d_{\textrm{def}}$
otherwise. The latter is computed as the average tether bond length in the conformation.

## Parameter file entry

```yaml
# Parameters for the tether forces.
tether forces:
  # Whether it's enabled.
  enabled: true
  H1: 100 eps/A**2
  H2: 0 eps/A**4
```