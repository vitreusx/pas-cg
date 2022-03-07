# Heurestic dihedral angle FF

The heurestic dihedral angle force field is described by:

$$ V = \sum f(\varphi_i, \mathtt{type}_i)\\ f(\varphi, \mathtt{type}) = a + b
\sin\varphi + c \cos\varphi + d\sin^2\varphi + e\cos^2\varphi + f\sin\varphi
\cos\varphi $$

where the summation is over the quadruples without the native value and the
coefficients are different for each $\mathtt{type}$. The type of the quadruple
is determined in the same fashion as with the heurestic dihedral types, with the
second and the third of the four residues determining the type.

## Parameter file entry

```yaml
heurestic dihedrals:
  enabled: boolean
  coefficients: csv file
```

The coefficients csv file has the
schema `type2,type3,const,sin,cos,sin2,cos2,sin_cos`, where the coefficients (1)
correspond to $a, \ldots, f$ as above, (2) are assumed to be in $\varepsilon$ if
given without a unit. The values for all 9 combinations of `type2` and `type3`
must be specified.