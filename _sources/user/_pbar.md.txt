# Progress bar

No particular explanation is necessary. The data displayed includes the progress
to `general.total time`, wall clock time and values of $V$ and $t$. The progress
bar is updated based on the wall clock time.

## Parameter file entry

```yaml
# Parameters for the progress bar.
progress bar:
  # Whether it's enabled.
  enabled: true

  # Width of the progress bar, in characters.
  width: 64

  # Real-time update period.
  update period: 2.5 s
```