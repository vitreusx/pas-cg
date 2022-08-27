# Native contacts FF

The native contacts force field is described by:

$$ V = \sum \mathsf{LJ}(r, r_0) + \sum \mathsf{Dis}(r, r_0)
$$

where the first summation is over the native contacts which are not disulfide
bonds, the second summation is over the disulfide bonds and $r_0$ is the native
distance between the residues in contact. $\mathsf{LJ}$ denotes the standard
Lennard-Jones potential:

$$ \mathsf{LJ}(r, r_0) = \varepsilon
\left[\left(\frac{r_0}{r}\right)^{12}-2\left(\frac{r_0}{r}\right)^6\right]
$$

```{note}
Note: The depth of the LJ potential is the same for each non-disulfide native contact, and equal to $1 \varepsilon$. The different types of the contacts (like backbone-backbone or sidechain-sidechain) are also not distinguished.
```

A special potential can be specified for the native disulfides (i.e. such as
appear in the `SSBOND` pdb file entries, **not** the derived cysteine pairs or
such as may be specified in the sequence file).

(foldingunfolding-study)=
## Folding/unfolding study

## Parameter file entry

```yaml
# Parameters for the Go-model contacts, as well as some "simulation protocol"
# settings pertaining to these.
native contacts:
  # Whether it's enabled.
  enabled: true

  # Depth of the LJ potential used.
  lj depth: 1 eps

  # Parameters for the native disulfide bonds (the SSBOND entries in
  # the PDB file). NOTE: I'm not sure if it should include all native contacts
  # two cysteines, or only the ones for which an SSBOND record exists.
  disulfide bond force:
    variant: harmonic
    harmonic params:
      H1: 100 eps/A**2
      H2: 0 eps/A**4
      nat_r: 6.0 A
    lj params:
      depth: 4.0 eps
      r_min: 6.0 A

  # Parameters specifying the folding/unfolding study/protocol.
  (un)folding study:
    # Enable early stopping when all contacts are formed or broken (more
    # generally, having changed the state from the native one).
    stop when all are formed/broken: false

    # Measure median first-change (formed -> broken or broken -> formed for
    # unfolding/folding studies respectively) time across multiple trajectories.
    measure median times: false
```