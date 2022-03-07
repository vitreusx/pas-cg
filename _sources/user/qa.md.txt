# Quasi-adiabatic (QA) FF

The quasi-adiabatic (QA) force field is non-Hamiltonian, and is based on the
idea of switchable contacts. The algorithm has three components in total:

1. For the non-local pairs not in contact, we check whether they satisfy the _
   formation criteria_, and if they do, they are added to the list of _candidate
   contacts_.
2. For the non-local pairs in contact, we evaluate the potential resulting from
   the contact, and check if the contact satisfies the _breaking criteria_. If
   the contact is completely desaturated, we *mark* them as removed.
3. **Sequentially**, we remove the marked contacts, and process the candidates (
   in the order of increasing distances) by checking whether the contact can
   actually be formed, and if so, adding it to the list of contacts.

The split of pre-checking the formation criteria in parallel and doing it again
sequentially later on a sifted list is done so as to ensure determinism of the
algorithm without the checking becoming a bottleneck, assuming the candidate
contacts are few in number. As will be shown later, the formation and breaking
of contacts affects the creation of new contacts.

```{note}
For a given simulation step, the contacts formed in that step do not participate
 in the evaluation of the potential in the same step. Rather, sequential
  processing is done after the asynchronous part of the simulation step,
   which contains the potential evaluation.
```

## Contact

When contact is formed or broken, it does not assume its full strength
immediately, but rather increases linearly from no strength to full strength in
a given time span.

A QA contact can be described by the following data:

- participating residues;
- contact _type_: either backbone-backbone (bb), backbone-sidechain (bs, sb), a
  pair-type-specific version of the sidechain-sidechain (ss) contact. Thus,
  there are $N^2+3$ different types of QA contacts, where $N=20$ is the number
  of amino acid types;
- status: whether a contact is saturating (i.e. increasing in or maintaining
  strength) or desaturating (i.e. decreasing in strength);
- reference time: contact creation time if the contact is saturating, or contact
  breaking time if the contact is desaturating. Used for computing the
  saturation.

## Formation criteria

A QA contact between two residues can be formed, if three criteria are
satisfied:

### Distance criterion

To each type of the contact type (as described above) is associated a distance
$r_{\textrm{min}}$. The distance between the residues must be lower than $r_
{\textrm{min}}(1+\mathtt{formation\_tolerance})$ for the criterion to be
satisfied.

### Geometric criterion

We approximate the directions of the backbone hydrogen bond or the sidechain $C_
\beta$ atom from the positions of three consecutive $C_\alpha$ atoms. More
specifically, for a consecutive triple, we may define:

$$ \vec{n_2} = \frac{\vec{r_{23}} - \vec{r_{12}}}{|\vec{r_{23}} - \vec{r_
{12}}|}\\ \vec{h_2} = \frac{\vec{r_{23}} \times \vec{r_{12}}}{|\vec{r_{23}}
\times \vec{r_{12}}|} $$

which are vectors associated with the central site. The direction of the
sidechain is approximately $-\vec{n}$, and the direction of the hydrogen bond is
$\pm \vec{h}$. Based on these, we establish following criteria for the different
types of contacts:

- bb: $|\cos(\vec{h_1}, \vec{r_{12}})|, |\cos(\vec{h_2}, \vec{r_{21}})| \geq
  \mathtt{min\_abs\_cos\_hr}$ and $|\cos(\vec{h_1}, \vec{h_2})| >
  \mathtt{min\_abs\_cos\_hh}$;
- bs: $|\cos(\vec{h_1}, \vec{r_{12}})| \geq \mathtt{min\_abs\_cos\_hr}$ and
  $\cos(\vec{n_2}, \vec{r_{21}}) < \mathtt{max\_cos\_nr}$;
- sb: analogous to the bs case;
- ss: $\cos(\vec{n_1}, \vec{r_{12}}), \cos(\vec{n_2}, \vec{r_{21}}) <
  \mathtt{max\_cos\_nr}$.

where $\cos(\vec{u}, \vec{v})$ denotes the cosine of the angle between the two
vectors. In general, $|\cos(\vec{h_i}, \vec{r_{ij}})| \geq
\mathtt{min\_abs\_cos\_hr}$ is the criterion for the backbone of residue $i$,
and $\cos(\vec{n_i}, \vec{r_{ij}}) < \mathtt{max\_cos\_nr}$ is the criterion for
the sidechain of residue $i$, the extra condition for the bb case
notwithstanding.

### Sync values criterion

With each residue are associated "sync values", which are numbers which
essentially track the remaining "slots" for the contacts, in effect limiting the
number of contacts that can be formed. There are four such numbers:

1. Backbone contact slots.
2. Total sidechain contact slots.
3. Sidechain contact slots for the _polar_ residues.
4. Sidechain contact slots for the _hydrophobic_ residues.

The base values are amino acid-specific. The presence of this criterion is the
cause for the split of pre-filtering the candidate contacts and processing them
sequentially later, as the formation of a contact could reduce the sync values
to such a state as to prevent the formation of another contact in another
thread, which would cause the algorithm to be nondeterministic. Similarly, we
cannot "just" remove the desaturated contacts, but must rather mark them first
and remove them sequentially later.

```{note}
The _native_ contacts do not affect the sync values.
```

## Force field

The residues in contact interact with a regular Lennard-Jones potential, unless
sinking version of the ss-type contacts is turned on (see LJ variants section
for more details), in which case it is used instead. This potential is then
multiplied by a _saturation factor_, which increases from $0$ to $1$ in time
span of $\mathtt{phase\_duration}$ when the contact is formed, and from $1$ to
$0$ when the contact is broken (in the same time span).

```{note}
Time derivative of the saturation is assumed to be zero in the computation.
```

## Breaking criterion

The only breaking criterion is for the distance between the residues to become
greater than $\mathtt{breaking\_factor} \cdot 2^{-1/6} r_{\textrm{min}}$.

## Dynamic disulfide bonds

The procedure for the pairs of cysteines can be altered so as to simulate the
formation of dynamic disulfide bonds. This is accomplied in two orthogonal ways:

1. The formation and breaking criteria can be altered.
    - For a disulfide bond to be created between two cysteines, the following
      criteria must be satisfied:
        1. The usual criteria for an ss-type contact must be satisfied.
        2. The cysteines in question have not already formed disulfide bonds
           with other cysteines
        3. The sum of the number of neighbors (i.e. number of residues within
           $\mathtt{neigh\_radius}$) is smaller than
           $\mathtt{max\_neigh\_count}$.
    - For a disulfide bond to be broken, following criteria must be satisfied:
        1. The distance $d$ is such, that $|d - \mathtt{def\_bond\_dist}| \geq
           \mathtt{max\_dist\_dev}$.
        2. The sum of the number of neighbors is smaller than
           $\mathtt{max\_neigh\_count}$.
2. The force can be altered, setting it either to be harmonic, or of a
   Lennard-Jones type.

```{warning}
At the moment, the neighbor counter only counts the "nonlocal" pairs.
```

## Parameter file format

```yaml
quasi-adiabatic:
  enabled: boolean
  phase duration: quantity [T]
  formation tolerance: number
  breaking factor: number
  min |cos(h, r)|: number
  min |cos(h, h)| for bb: number
  max cos(n, r): number
  disulfides:
    force: disulfide force
    special criteria:
      enabled: boolean
      max neigh count: integer
      neigh radius: quantity [L]
      def bond dist: quantity [L]
      max bond dist deviation: quantity [L]
```

The `disulfide force` type of entry is as in the "Native contacts" section.

## Output data

Following data is emitted per snapshot:

- A list of current sync values, into `sync_values.csv`;
- A list of current contacts, into `qa_contacts.csv`. Aside from the
  aforementioned data per contact, computed saturation is also included.
- Aggregate stats, similar to the ones for the native contacts. Aside from the
  ones enumerated there, we also include the numbers for the dynamic disulfide
  bonds. 