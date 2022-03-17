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

## Parameter file entry

```yaml
native contacts:
  enabled: boolean
  lj depth: quantity [E]
  active threshold: number
  disulfide bond force: disulfide force

disulfide force:
  harmonic:
    H1: quantity [E/L**2]
    H2: quantity [E/L**4]
    equilibrium dist: quantity [L]
  lj:
    depth: quantity [E]
    r_min: quantity [L]
```

The `disulfide bond force` entry is optional - if present, the special potential
is enabled. The `harmonic` and `lj` entries of the `disulfide force` entry are
mutually exclusive.

## Output data

Following data is emitted per snapshot:

- a list of all native contacts, into `nat_conts.csv`. The columns of the file
  are:
    - `i1`, `i2`: the indices of the residues in contact;
    - `chain1`, `chain2`: corresponding chain indices, for the purposes of
      determining whether the contact is intra- or interchain;
    - `nat_dist[A]`, `cur_dist[A]`: native and current distance (in Angstrem);
    - `is active`: whether a contact is currently active;
    - `formed once`: whether the contact has been formed previously;
    - `formation time`: if `formed once` is `true`, the (internal) time when it
      was formed for the first time.
- a number of "aggregate" statistics derived from this list for the active
  contacts, placed into the `num of active contacts` node:
    - `all`: number of all active contacts;
    - `back-back`: number of backbone-backbone contacts (both inter- and
      intra-chain);
    - `back-back (intra)`, `back-back (inter)`: number of bb contacts with given
      type;
    - analogous statistics for `back-side` and `side-side` type, also with the
      intra- and inter-chain distinction.

Some of these statistics are added to the `scalars.csv` file, namely:

- `nc-all` for all contacts;
- `nc-b1-b1` for intra-chain backbone-backbone contacts;
- `nc-b1-b2` for inter-chain backbone-backbone contacts;
- `nc-b1-s1`, `nc-b1-s2`, `nc-s1-s1`, `nc-s1-s2`, analogous.

```{warning}
The data is emitted unconditionally whenever the potential is enabled.
```