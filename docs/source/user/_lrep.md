# Local repulsive force field

The local repulsive FF acts on residues $(i, i+2)$ in the same fashion as [the Pauli exclusion potential](_pauli).

## Parameter file entry

```yaml
# Parameters for local repulsive potential.
local repulsive:
  # Whether it's enabled.
  enabled: true

  # Depth of the potential.
  depth: 1.0 eps
```
