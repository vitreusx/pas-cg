# Chirality FF

The chirality force field is described by:

$$ V = \sum \frac{1}{2} e_{\mathrm{chi}} (C-C_0)^2\\ C_{1234} = \frac{\langle
\vec{r_{12}} \times \vec{r_{23}}, \vec{r_{34}} \rangle}{r_{0,23}^3} $$

where $r_{0,ij}$ is $r_{ij}$ in the native structure, and $C_0$ is the value of
$C$ in the native structure. The summation is over all the quadruples. The force
field requires the model to have been loaded from the PDB file, since the
sequence file lacks the values of the native tether distances.

## Parameter file entry

```yaml
chirality:
  enabled: boolean
  e_chi: quantity [E]
```