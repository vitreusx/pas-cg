# Amino acid data

Various properties of the amino acid types can be specified in the form of a
YAML file.

```yaml
amino acid data:
  default atom data:
    (atom name): quantity [L]

  amino acids:
    (amino acid name):
      mass: quantity [M]
      alt atom radii:
        (atom name): quantity [L]
      side: list of atom names
      radius: quantity [L]
      polarization: one of [hydrophobic, polar]
      charge: quantity [(dim of e)]
      contact limits:
        back: integer
        side (all): integer
        side (hydrophobic): integer
        side (polar): integer
```

The fields `alt atom radii`, `polarization` and `charge` are optional. For more
info, see the `amino acids` entry in
the [`default.yml`](https://github.com/vitreusx/pas-cg/blob/main/cg/src/simul/default.yml)
file.