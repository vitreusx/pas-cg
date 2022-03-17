# Program output

The output of the program consists solely of the progress bar and the actual
data, emitted into files placed in the output directory. The data is emitted at
regular intervals (in terms of internal time).

The structure of the output directory is as follows:

```
[output]
├── report.yml
├── [other per-simul files...]
├── traj-${traj_idx}
|   ├── report.yml
|   ├── scalars.csv
|   ├── model.pdb
|   ├── [other per-trajectory files...]
|   ├── snap-${snap_idx}
|   |   ├── report.yml
|   |   ├── [other per-snapshot files...]
|   |   |   ...
```

The data is outputted to the `report.yaml` files, as well as other files (for
example, `model.pdb` would contain the PDB model of the system). The data is
organized into three levels:

- data pertaining to the entire simulation (i.e. all trajectories);
- data pertaining to a given trajectory;
- data pertaining to a given *snapshot*, emitted at a regular interval. This
  report may also include files - however, there is an option to emit scalar
  statistics and the entire files at different intervals.

```{warning}
At the moment, one does **not** have the option for including or excluding data from the
output, so for example one may not choose to skip computing the RMSD.
```

## Per-simulation report

At the moment, no per-simulation information is reported.

```{note}
One could output the `inputfile.yml` here, in a fashion altogether similar to the original headers of the .out file.
```

## Per-trajectory report

For a given trajectory, following data is emitted:

- `model.pdb` file, containing `MODEL` entries for all the snapshots. This is
  done chiefly in order to be able to quickly investigate the evolution of the
  system in a tool such as [PyMOL](https://pymol.org/2/).
- `scalars.csv` file, containing a table with selected scalar values from all
  the snapshots.

## Per-snapshot report

For a given snapshot, following data is emitted:

- "general" statistics: time, potential, kinetic and total energy;
- positions of the monomers in the form of a PDB file, `model.pdb`;
- current values of tethers, bond and dihedral angles, along with native ones if
  they exist, placed into the files `tethers.csv`, `angles.csv`
  and `dihedrals.csv`;
- "kinematic" properties: $\lambda_x$, $\lambda_y$, $\lambda_z$, radius of
  gyration, asphericity (see the article on
  the [gyration tensor](https://en.wikipedia.org/wiki/Gyration_tensor) for more
  details);
- `chains.csv` file with the kinematic properties for every chain separately;
- data pertaining to [native contacts](nat_cont.md);
- data pertaining to [quasi-adiabatic potential](qa.md) and the created
  contacts;
- data pertaining to the [AFM tips](afm.md).

The scalars added to the `scalars.csv` file are:

- `t`: time;
- `V`: potential energy;
- `K`: kinetic energy;
- `E`: total energy;
- `Rg`: radius of gyration of the entire structure;
- `W`: asphericity of the entire structure.

```{warning}
Asphericity is, from what I've seen, computed in a slightly different fashion 
than what's in the wiki (it seems to be scaled in some sense). This particular
implementation follows the Wiki article.
```

## Parameter file entry

```yaml
output:
  enabled: boolean
  stats period: quantity [T]
  file period: quantity [T]
```
