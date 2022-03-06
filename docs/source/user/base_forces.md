# Base forces. Variants of the Lennard-Jones potential

Different potentials can be set for backbone-backbone, backbone-sidechain and
sidechain-sidechain contacts. For the first two types, the only option is the
standard Lennard-Jones potential, with adjustable depth and $r_{\textrm{min}}$.
For sidechain-sidechain contacts, there are more options. First of all, the
parameters can be either set by default, or adjusted for every combination of
residue types. Second, a _sinking_ Lennard-Jones potential can be used.

$$ \mathsf{SinkLJ}(r) = \begin{cases} \mathsf{LJ}(r) & r \leq r_{\textrm{min}}\\
\mathsf{LJ}(r_{\textrm{min}}) & r_{\textrm{min}} \leq r < r_{\textrm{max}}\\
\mathsf{LJ}(r_{\textrm{min}} + r - r_{\textrm{max}}) & r_{\textrm{max}} \leq r
\end{cases} $$

In other words, it's a normal Lennard-Jones potential, but where the potential
well has been expanded from just $r_{\textrm{min}}$ so as to span an interval
$[r_{\textrm{min}}, r_{\textrm{max}}]$.

## Parameter file entry

```yaml
lj force variants:
  bb:
    r_min: quantity [L]
    depth: quantity [E]
  bs:
    r_min: quantity [L]
    depth: quantity [E]
  ss:
    use sinking variant: boolean
    default:
      r_min: quantity [L]
      r_max: quantity [L]
      depth: quantity [E]
    per-pair:
      r_min: csv file
      r_max: csv file
      depth: csv file
```

The different options in `ss.default` and `ss.per-pair` can be omitted or
included, though for every parameter one of them should be included.