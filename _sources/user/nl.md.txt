# Neighbor list

In the current version, the neighbor list is in fact a sorted list of pairs of
residues within a given cutoff distance. This cutoff distance is the maximum of
the cutoff distances for the non-local force fields which use the neighbor list.
In order to not recompute the list every time, a pad is introduced - when (re)
computed, we include all pairs of residues within the distance of
$\mathtt{cutoff}+\mathtt{pad}$ of each other, and at every step we verify
whether **the sum** of:

1. twice the maximal displacement of the residues from their positions when the
   list was recomputed;
2. $L^1$ displacement of the simulation box size;

exceeds the value of the pad.

Two algorithms for the computation exist:

1. The "legacy" version, which goes over all the pairs and checks the distance;
2. The "cell" version, which first assigns the residues to the cells, and goes
   over the pairs of residues in the same or neighboring cells.

The algorithm computes two different lists: of the non-native pairs, and of the
native pairs.

## Parameter file entry

```yaml
neighbor list:
  algorithm: one of [legacy, cell]
  pad: quantity [A]
```