# Input/parameter file

The execution of the program is controlled by modifying the input file provided in the `-p/--param-file` option.

(default_yml)=
## `default.yml` - Default input file

The reference input file with all the possible parameters and their default values, along with comments, is `cg/default.yml` - if in doubt, you can consult this file. You can also use this file as a template, though it is not necessary - internally, the two files are merged, with the user input file overriding the values in the default, so that, for example, user input file can be empty, or contain only those values which are changed. You can see the example input files in `cg/data` for more details. It is, moreover, highly advisable to read the default file whilst consulting the documentation on the various available settings, sections etc.

## Structure of the input file

The input file is in a YAML format, extended by the possibility of referencing files - see [](_formats) for further details. All the parameters are organized into sections, also called stanzas.

### Force specification

In a number of places, one inserts values for a given type of a "base" force field, like LJ, shifted LJ, harmonic etc. This format is as follows:

```yaml
# An example spec for the regular forces.
a random force:
  # Force variant to be used. Options:
  # - "lj": Lennard-Jones force;
  # - "sink lj": a "sinking" LJ force, in which the pointlike minimum of the
  #   potential is extended to span a certain interval of the radii.
  # - "shifted lj": a version of the LJ force, effective only at close distances
  #   (i.e. below r_min) and shifted to make the value of the potential always
  #   positive.
  # - "harmonic": harmonic force.
  variant: lj

  # Parameters for the "lj" option.
  lj params:
    # Depth of the potential.
    depth: 1 eps

    # Distance value, at which the potential reaches minimum; equivalently,
    # 2^{1/6} \sigma.
    r_min: 5 A

  # Parameters for the "shifted lj" option - the same, as for the regular LJ
  # potential.
  shifted lj params:
    depth: 1 eps
    r_min: 4 A

  # Parameters for the "sink lj" option.
  sink lj params:
    # Depth of the potential, as in the regular case.
    depth: 1 eps

    # The start of the "sink" of the potential, i.e. the interval of distances
    # at which the potential reaches its minimum.
    r_low: 5 A

    # The end of the aforementioned interval.
    r_high: 9 A

  # Parameters for the "harmonic" option.
  harmonic params:
    # Value of H1 in V = H1 d^2 + H2 d^4
    H1: 100 eps/A**2

    # Value of H2 in V = H1 d^2 + H2 d^4
    H2: 0 eps/A**4

    # Optional "blanket" value of the native distance, used for example with
    # disulfide bonds.
    nat_r: 6.0 A

a sidechain-sidechain force:
  variant: lj
  lj params:
    # For each force parameter, we can set the "default"
    # per-residue-pair parameters.
    default:
      depth: 1 eps

    # Similarly, we can provide a table of values, in a format similar to the
    # table of minima of the potentials in `ss_pairs_csv` entry.
    per pair:
      r_min:
        (at path): values.csv
```

See [](default_yml) for more details.

### Internal units and quantities.

All the physical quantities in the input data can be given units. The quantity
format is either `[number]`, in which case it is said to be dimensionless,
or `[numerical value] [unit]`.

The list of predefined units is: "f77unit" (internal unit of length, equal to 5
Angstrem), "A" (angstrom), "nm", "m", "ns", "tau" (internal unit of time, equal
to 1 ns), "micros" (microseconds), "ms", "s", "atom", "mol", "eps" (internal
unit of energy, equal to $\approx 1.5 \textrm{kcal}/\textrm{mol}$), "kcal", "J"
, "kB" (Boltzmann constant), "K", "kg", "amu", "f77mass" (internal unit of mass,
equal to average residue mass), "e", "C", "Amp" (Ampere, since "A" is taken by
angstrom), "c", "H", "mu_0" (magnetic permittivity of the vacuum), "eps_0" (
electric permittivity of the vacuum), "rad", "deg".

The units can be divided or raised to nth power (with `**`) - in general, any
arithmetic expression involving these predefined units can be entered.

```{warning}
The system of units is inconsistent. From `f77unit`, `tau` and `eps` one could derive the base unit of mass, yet `f77mass` is **not** equal to it.
```

For more info,
see [`units.h`](https://github.com/vitreusx/pas-cg/blob/main/cg/include/cg/utils/units.h)
and [`quantity.cpp`](https://github.com/vitreusx/pas-cg/blob/main/cg/src/utils/quantity.cpp)
.

## Stanzas

Below, you can find the description of all the sections/stanzas available. Note that for many of them, we attach links to pages dedicated to the given functionality.

### `general` - General settings

This stanza contains settings which do not belong to any more specialized section.

```yaml
# General settings
general:
  # Program mode. Options:
  # - "perform simulation": Standard model
  # - "check determinism": Runs two copies of simulation in parallel to see
  #   if the results are different. Effective when # of threads > 1; one can
  #   check in this way if the results obtained with multi-thread simulations
  #   are reproducible.
  mode: perform simulation

  # Total, maximal time of simulation for a single trajectory. The actual
  # execution time can be lower; in the other direction, processes such as
  # oscillations can be prematurely terminated in this fashion.
  total time: 3e6 tau

  # Equilibration time, i.e. time spent before "external mechanisms", such as
  # pulling chains by AFM residues or squeezing, are turned on. The forces of
  # friction and the thermostat are still active during equilibration, though.
  equil time: 0 tau

  # Simulation seed.
  seed: 448

  # Number of threads used.
  num of threads: 1

  # Number of trajectories sampled.
  num of trajectories: 1

  # Cutoff for the repulsive interactions. This includes local repulsive forces
  # and Pauli exclusion forces.
  repulsive cutoff: 5 A

  # Counting factor, i.e. multiplier for \sigma in LJ potential, below which
  # the contact is considered to be "active".
  counting factor: 1.3

  # Debug options
  debug mode:
    # Turn floating point exceptions on/off.
    floating point exceptions: false
    # Dump raw data at every step - used in the parity tests.
    dump data for every step: false
```

### `input` - [](_input)

### `simulation box` - [](_sbox)

### `angle potentials` - [](_angle)

### `Debye-Hueckel` - [](_dh)

### `pseudo-improper dihedral` - [](_pid)

### `quasi-adiabatic` - [](_qa)

### `chirality` - [](_chir)

### `native contacts` - [](_nc)

### `Pauli exclusion` - [](_pauli)

### `tether forces` - [](_tether)

### `progress bar` - [](_pbar)

### `output` - [](_output)

### `checkpoints` - [](_ckpt)

### `local repulsive` - [](_lrep)

### `amino acid data` - [](_amino)