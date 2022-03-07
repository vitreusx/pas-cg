# Atomic Force Microscope (AFM)

We can simulate the actions of pulling and stretching of the proteins by an
Atomic Force Microscope (AFM). There are two options available regarding the
type of the AFM:

- pulling with constant velocity, starting from the original position of the
  pulled residue and moving in some direction with a given velocity;
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
afm:
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