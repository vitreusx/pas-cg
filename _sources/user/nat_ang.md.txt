# Native bond angle FF

The native bond angle force field is described by:

$$V = \sum \frac{1}{2}k(\theta - \theta_0)^2$$

where the summation is over the triples with the native value of the bond angle
present.

The bond angle itself is determined as follows:

$$\theta_{123} = \arccos
\left[-\frac{\langle \vec{r_{12}}, \vec{r_{23}} \rangle}{r_{12} r_{23}}\right]$$

## Parameter file entry

```yaml
native angles:
  enabled: boolean
  k: quantity [E/Angle^2]
```