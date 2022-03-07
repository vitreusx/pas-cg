# Heurestic bond angle FF

The heurestic bond angle force field is described by:

$$ V = \sum f(\theta_i, \mathtt{type}_i)\\ f(\theta, \mathtt{type}) = \sum_
{i=0}^D c_i[\mathtt{type}] \theta^D $$

where the summation is over the triples with the native value of the bond angle
not present. In other words, we approximate the value of the potential by a set
of polynomials, one for each different _type_ of the triple of residues. The
type of the triple is derived in the following fashion:

1. First, we class the residues into glycines (code `G`), prolines (code `P`)
   and all other amino acids into a single category (code `X`);
2. Then, the type of the triple of residues corresponds to the pair consisting
   of the class of the second and of the third residue. Thus for example, if we
   have a triple of residues which are glycine, proline and valine, the type of
   the residue is `PX` (`P` for the proline, `X` for the valine).
3. Thus, there are nine total types: `GG`, `GP`, `GX`, ..., `XP`, `XX`.

## Parameter file entry

```yaml
heurestic angles:
  enabled: boolean
  coefficients: csv file
```

The coefficients csv file has the schema `type1,type2,a0,a1,a2,a3,a4,a5,a6`,
where the coefficients are assumed to be in $\varepsilon/\mathrm{rad}^d$ for
$\mathtt{a}d$ if given without a unit. The file should contain entries for all 9
combinations of `type1` and `type2`.

```{warning}
The interpolation is such that, for the "non-regular" values of $\theta$, such 
as may occur in the beginning of the simulation after morphing into a 
self-avoiding walk, the value of $f(\theta, \mathtt{type})$ may be extremely large.
In order to prevent any issues, the magnitude of the force is clamped to 
$[-10^3, 10^3]$ (in internal units).
```