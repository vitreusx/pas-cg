# General parameters

```yaml
general:
  total time: quantity [T]
  equil time: quantity [T]
  seed: integer
  num of threads: integer
  debug mode:
    enabled: boolean
    floating point exceptions: boolean
    determinism: boolean
  num of trajectories: integer
  disable all forces: boolean
```

For any of the settings in the debug mode, `general.debug mode.enabled` must be
set to `true`.

`general.debug mode.determinism` flag turns on a special mode of the program, in
which two programs are run in parallel and compared against each other. For this
mode, output is disabled.

`general.disable all forces`, well, disables all forces, except for the Langevin
thermal noise.