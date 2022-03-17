# Atomic Force Microscope (AFM)

We can simulate the actions of pulling and stretching of the proteins by an
Atomic Force Microscope (AFM). There are two options available regarding the
type of the AFM:

- pulling with constant velocity, starting from the original position of the
  pulled residue and moving in some direction with a given velocity

```{note}
The *AFM tip* has the constant velocity, not the residue itself.
```

- pulling with constant force, which (in effect) adds a given force vector to
  $F_i$.

As for specifying the pulled residues, we may choose to either:

- pull a single residue in a given direction;
- pull terminal residues of a chain apart (in the direction determined by the
  original positions of the residues in question) with a force or velocity of a
  given magnitude (since the direction is determined).

The tip of the AFM is connected with the pulled residue by a harmonic tether
with given values of $\mathtt{H1}$ and $\mathtt{H2}$.

## Parameter file entry

```yaml
AFM:
  H1: quantity [E/L**2]
  H2: quantity [E/L**4]
  tips:
    - type: one of [const velocity, const force]
      single residue: integer
      direction: [ x, y, z ]
    - type: (as above)
      pulled-apart chain: integer
      magnitude: quantity (of force or of velocity, depending on the type)
```

For the "single residue" specification, the direction should be a vector of
quantities appropriate for the type of the tip.

## Output

Three files are emitted regarding the state of the AFM tips:

- `vel-afm.csv`: state of the const-velocity AFM tips. For a pulled-apart
  chains, two entries (one for each terminal residue) is included. The columns
  are:
    - `res idx`: index of the pulled residue;
    - `force [f77unit*A/tau**2]`: force acting on the residue.
- `force-afm.csv`: analogous to `vel-afm.csv` but for the const-force AFM tips.
  The columns are:
    - `res idx`: as above;
    - `vel [A/tau]`: velocity of the residue.
- `pulled-chains.csv`: extra data for the pulled-apart chains:
    - `chain idx`: the index of the chain;
    - `end-to-end distance [A]`: the distance between the terminal residues.

```{warning}
At the moment, the forces and the velocities are **not** averaged.
```

If there is only one pulled chain, a scalar is added to
the `scalars.csv`: `chain length [A]`, i.e. the `end-to-end distance [A]` for
the sole pulled-apart chain.