# Program output

The output of the program consists solely of the progress bar and the actual
data, emitted into files placed in the output directory. The data is emitted at
regular intervals (in terms of internal time).

The structure of the output directory is as follows:

```
[output]
├── report.yml
├── [other files...]
├── [snapshot number]
    ├── report.yml
    ├── [other files...]
```

In `report.yml`, one can find data pertaining to the simulation as a whole.
Beyond that, data pertaining to the state of the simulation in a given moment is
emitted in the form of zero-indexed snapshots.

## Emitted data and files

For the simulation as a whole, following data is emitted:

- `model.pdb` file, containing `MODEL` entries for all the snapshots. This is
  done chiefly in order to be able to quickly investigate the evolution of the
  system in a tool such as [PyMOL](https://pymol.org/2/).

At regular intervals, following data pertaining to the state of the simulation
in a given moment is emitted:

- "general" statistics: time, potential, kinetic and total energy;
- positions of the monomers in the form of a PDB file, `model.pdb`;
- current values of tethers, bond and dihedral angles, along with native ones if
  they exist, placed into the files `tethers.csv`, `angles.csv`
  and `dihedrals.csv`;
- "kinematic" properties: $\lambda_x$, $\lambda_y$, $\lambda_z$, radius of
  gyration, asphericity (see the article on
  the [gyration tensor](https://en.wikipedia.org/wiki/Gyration_tensor) for more
  details);
- data pertaining to native contacts;
- data pertaining to quasi-adiabatic potential and the created contacts.

```{warning}
Asphericity is, from what I've seen, computed in a slightly different fashion 
than what's in the wiki (it seems to be scaled in some sense). This particular
implementation follows the Wiki article.
```

## Parameter file entry

```yaml
output:
  enabled: boolean
  period: quantity [T]
```

```{warning}
One does **not** have the option for including or excluding data from the
output, so for example one may not choose to skip computing the RMSD.
```