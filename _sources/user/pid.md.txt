# Pseudo-improper dihedral (PID) force field

Pseudo-improper dihedral (PID) force field is the second custom potential. It
acts on "bundles" which consist of two central, interacting, residues (
numbered `2` and `5`), along with the neighboring residues (`1` and `3` for the
residue `2` and `4`, `6` for the residue `5`). The general form of the potential
is:

$$ V = \sum V_b(\mathtt{bundle})
$$

where the sum is over all the bundles, i.e. all pairs of non-terminal
non-natively-connected residues along with the neighbors. The bundle term is:

$$ V_b = \sum_i \lambda_i(\psi_{24}) \lambda_i(\psi_{42}) \phi_i(r_{24})
$$

The sum runs over three components of the PID potential: "bb+", "bb-" and "ss"
components. These differ in the choice of the Lennard-Jones potential $\phi_i$,
as well as the parameters of the $\lambda$ functions.

The $\lambda$ functions are Gaussian-like, and are parametrized by two
parameters: $\alpha$ and $\psi_0$. There are two variants of the mode of
evaluation of the $\lambda$ functions:

1. The "cosine" $\lambda$ function:

   $$ \lambda(\psi) = \begin{cases} \frac{1}{2}(1+\cos \alpha(\psi-\psi_0)) &
   |\alpha (\psi - \psi_0)| < \pi\\ 0 & \textrm{otherwise} \end{cases} $$

2. The "algebraic" $\lambda$ function:

   $$ \lambda(\psi) = \frac{t^2-2|t|+1}{2t^2-2|t|+1}, t = \frac{\alpha (\psi -
   \psi_0)}{\pi} $$

The terms $\psi_{42}$ and $\psi_{24}$ are the improper dihedral angles. For a
configuration of three consecutive residues `123` and another residue `4`, by
$\psi_{24}$ we shall denote the angle between the planes `123` and `143`.

## Parameter file entry

```yaml
pseudo-improper dihedral:
  enabled: boolean
  lambda version: one of [cosine, algebraic]
  bb+:
    alpha: quantity [1/Angle]
    psi_0: quantity [Angle]
    r_min: quantity [L]
    depth: quantity [E]
  bb-: as for bb+
  ss:
    alpha: quantity [1/Angle]
    psi_0: quantity [Angle]
```