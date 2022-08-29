# Pseudo-improper dihedral (PID) force field

Pseudo-improper dihedral (PID) force field is the second custom potential. It
acts on "bundles" which consist of two central, interacting, residues (
numbered `2` and `5`), along with the neighboring residues (`1` and `3` for the
residue `2` and `4`, `6` for the residue `5`). The general form of the potential
is:

$$ V = \sum V_b(\mathtt{bundle})
$$

where the sum is over all the bundles, i.e. all pairs of non-terminal
non-natively-connected residues along with the neighbors; whether the residues
separated by 3 others should be included is controlled with
the `include separated by 3` setting. The bundle term is:

$$ V_b = \sum_i \lambda_i(\psi_{24}) \lambda_i(\psi_{42}) \phi_i(r_{24})
$$

The sum runs over two components of the PID potential: "bb+/-", and "ss"
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
# Parameters for the pseudo-improper dihedral (PID) potential.
pseudo-improper dihedral:
  # Whether it's enabled.
  enabled: false

  # Whether to allow for interactions between residues separated by three other
  # residues.
  include (i, i+4): true

  # Parameters for the "lambda" functions.
  lambda:
    # Variant of the lambda function used: Options: "algebraic", "cosine"
    variant: algebraic

    # Function parameters for the "bb plus" term.
    bb+:
      alpha: 6.4 1/rad
      psi_0: 1.05 rad

    # Function parameters for the "bb minus" term.
    bb-:
      alpha: 6.0 1/rad
      psi_0: -1.44 rad

    # Function parameters for the "ss" term.
    ss:
      alpha: 1.2 1/rad
      psi_0: -0.23 rad

  # Parameters for the forces.
  forces:
    # Variant of the forces used; in this case, values other than "lj" and
    # "sink lj" aren't supported. The case for the "sink lj" is somewhat
    # special, as by default it *inherits* the values from the lj params spec
    # below to r_high, and sets r_low to zero. In other words, if sinking lj
    # is enabled, by default the sink extends to zero distance, and the end
    # of the sink is as set in `ss_pairs_csv`.
    variant: lj

    # Parameters for the "bb plus" LJ potential.
    bb+:
      lj params:
        depth: 0.2 eps
        r_min: 5.6 A

    # Parameters for the "bb minus" LJ potential.
    bb-:
      lj params:
        depth: 0.2 eps
        r_min: 6.2 A

    # Parameters for the "ss" LJ potential.
    ss:
      lj params:
        default:
          depth: 1.0 eps
        per pair:
          r_min: *ss_pairs_csv
```