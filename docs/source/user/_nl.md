# Neighbor list

In the current version, the neighbor list is in fact a sorted list of pairs of residues within a given cutoff distance. In order to not recompute the list every time, a pad is introduced - when (re)computed, we include all pairs of residues within the distance of $\mathtt{cutoff}+\mathtt{pad}$ of each other, and at every step we verify whether **the sum** of:

1. twice the maximal displacement of the residues from their positions when the
   list was recomputed;
2. $\ell^1$ displacement of the simulation box size;

exceeds the value of the pad.

The algorithm computes two different lists: of the non-native pairs, and of the native pairs.

## Parameter file entry

```yaml
# Parameters for the neighbor (Verlet) list.
neighbor list:
  # Extra padding used to not have to recompute the list on each step.
  pad: 10.0 A

  # Overall cutoff of the interactions, excluding the electrostatics.
  cutoff: 20.0 A
```