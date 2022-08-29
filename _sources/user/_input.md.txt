# Input specifications

The program operates on protein chains, where each amino acid is represented by
its Ca atom. These chains can be provided to the program in one of two ways -
from a PDB file, or via a file written in custom sequence file format. Once that
is done, one can further modify the protein chains by arranging them into a
self-avoiding walk.

## Model data

The entire model, conceptually, consists of the following components:

1. Residues, with its amino acid type and an optional native position;
2. Chains, which are essentially lists of residues;
3. Native contacts, with the participating residues and the contact type. We
   distinguish following contact types: backbone-backbone, backbone-sidechain,
   sidechain-sidechain and disulfide bonds;
4. Pairs of connected residues, auto-generated from the chain structure. If the
   model was generated from the PDB file, an original length of the tether is
   saved;
5. Triples of connected residues, auto-generated from the chain structure. Each
   such triple angle may contain the native value of the bond angle, if such
   exists, in which case it is treated by the native bond angle force field.
   Otherwise, it is treated by the heurestic bond angle force field;
6. Quadruples of connected residues, auto-generated from the chain structure.
   Each such quadruple may contain the native value of the dihedral angle, if
   such exists, in which case it is treated by the complex or simple variant of
   the native dihedral angle force field. Otherwise, it is treated by the
   heurestic dihedral angle force field.

## Loading from a PDB file

A PDB file can be specified in the parameter file via either a path or an inline
source of the file. The protein chains in the file are reduced to the
coarse-grained representation by replacing the residues with the $C_\alpha$
atoms. Before that is done, one can derive a contact map between the residues
based on the distances between the $C_\alpha$ atoms or between all the atoms.
The contact is said to be established, if the distance between the atoms or
residues is less than the sum of radii of the atoms or residues in question
times $(26/7)^{1/6}$. Moreover, the native tether distances, bond angles and
dihedral angles are computed for the purposes of the relevant force fields.

## Loading from sequence file

A custom sequence file can be specified. Such a file consists of a list of
chains given by the amino acid codes of the consecutive residues. The chains are
by default intrinsically disordered, and therefore handled by the heurestic
force fields for the bond and dihedral angles. This can be overriden for
contiguous parts of the protein chains with the use of the map files. Such a map
file consists of the values of the native bond and dihedral angles, as well as a
list of native contacts between the residues.

## Sequence file format

Sequence file is written in a YAML format with a following "schema":

```yaml
model:
  chains: list of chain entries
```

where:

```yaml
chain:
  seq: string of amino acid codes for the consecutive residues
  maps: list of map entries

map: map-file subnode

map:
  __source: map-file subnode
  shift: integer

map-file:
  contacts: csv file
  angles: csv file
  dihedrals: csv file
```

The two formats of `map` entry (an unshifted one and a shifted one) are mutually
exclusive. The shift is added to all the indices present in the `map-file`.
The `map` entry of the `chain` node is optional - if missing, the chain is
intrinsically disordered.

Regarding the schemas of csv files:

- the schema of the contacts csv file is `i1,i2,length`, where `length` is
  assumed to be in Angstroms, if given without a unit;
- the schema of the angles csv file is `i1,i2,i3,theta`, where `theta` is
  assumed to be in radians, if given without a unit;
- thhe schema of the dihedrals csv file is `i1,i2,i3,i4,phi`, where `phi` is
  assumed to be in radians, if left without a unit.

## Morphing the chains into a self-avoiding walk

Positions of the residues can be set (or overriden in the case of the PDB file)
into a self-avoiding walk. The parameters for the procedure are:

- The initial residue density, from which an initial simulation box size is
  computed;
- Bond distance between the residues;
- Distance threshold for the residues to intersect;
- Number of retries;
- Whether the computed simulation box should be used in the simulation or
  discarded.

The procedure itself is as follows:

1. For each chain:
    1. We sample a position from the initial simulation box for the first
       residue in the chain, and sample a random direction in which the next
       residue will be placed.
    2. For the remaining residues, we place them at a bond distance from a
       previous residue in the direction of the current direction vector, and
       after that disturb the current direction vector by tilting it by an angle
       drawn uniformly from [$0$, $\pi/3$] and rotating it round the original
       value by an angle drawn uniformly from [$0$, $2\pi$]
2. We repeat it for all the chains, and after that check if any of the residues
   intersect:
    1. If they do, we start all over, for at most the given number of retries;
    2. Otherwise, we end the procedure;
    3. If the number of retries is exceeded, an exception is thrown.

## Parameter file stanza

```yaml
# Input specification
input:
  # Following two stanzas, "pdb file" and "seq file", are mutually exclusive,
  # and one of them is required for the program to function.

  # PDB file specification.
  pdb file:
    # Source file (verbatim or via path to file)
    #    source: ~
    #      (file contents)
    #      (at path): path to PDB file

    # Whether to load CRYST1 record from the PDB file, in which case the
    # simulation box size gets set to the size in the record.
    ignore CRYST1: false

    # Whether to skip unknown atoms in the PDB file.
    skip unknown: true

    # A map of aliases of the residues - some PDB files use nonstandard residue
    # names which are nominally not supported by the program.
    aliases: { HSD: HIS }

    # Sequence file specification
    seq file:
    # Source file (verbatim or via path to file)
  #      source: ~
  #        (file contents)
  #        (at path): path to seq file

  # Morphing the input model into a self-avoiding walk.
  morph into SAW:
    # Whether to perform the procedure
    perform: true

    # Bond distance used. In actuality, it's *slightly* more complicated, in
    # that this bond distance is used if the input comes from the sequence file,
    # and if the input comes from the PDB file, then the average of native
    # state tether lengths is used.
    bond distance: 3.8f A

    # A start box, i.e. the box from which the first residues in the chains are
    # sampled. The options are:
    # - scalar "origin", in which case the first residues are always given the
    #   position (0, 0, 0) - note that this **will** cause issues with more
    #   than one chain;
    # - a map with a "density" key, in which case the box spans from [-a/2,
    #   -a/2, -a/2] to [a/2, a/2, a/2], where a is computed so that the (molar)
    #   density of the start box is equal to the given value;
    # - a map with a "size" key, similar to the previous case but with an
    #   explicitly set box side length.
    start box: origin
    #    start box:
    #      density: 1e-4 residue/A**3
    #      size: 0 A

    # A distance, at which a given generated conformation is considered to be
    # self-intersecting, in which case we sample another conformation until
    # there are no such intersections.
    intersection at: 4.56 A

    # A number of retries of the procedure (with the retries being caused by
    # the aforementioned intersections).
    num of retries: 9001

    # Whether the procedure should take place in a start box with periodic
    # boundary conditions.
    periodic boundary conditions: false

  # Whether to use actual masses of the residues, as specified in the `amino
  # acid` stanza further ahead, or to average them out to a single value.
  normalize mass: true

  # Whether to load native structure (angles and dihedrals) from the PDB or
  # sequence file.
  load native structure: true
```