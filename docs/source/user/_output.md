# Output files

## Parameter file entry

```yaml
# Parameters for the output part of the program.
output:
  # Whether it's enabled.
  enabled: true

  # Write stats (.out file) every x time (in simulation time).
  emit stats every: 1e2 tau

  # Write structure (.pdb and .map files) every x time (in simulation time).
  emit structure every: 1e3 tau

  # Prefix for the output files (so, the files will be output.out, output.pdb
  # etc.)
  file prefix: output
```