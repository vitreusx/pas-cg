# Checkpoints

Checkpoints represent a complete state of the simulation at a given point. Indeed, simulation continued from a point of the creation of a checkpoint will proceed in exactly the same fashion as the original simulation in which the checkpoint was created.

## Loading the checkpoint

One can start the program from a checkpoint, using a `-c/--ckpt-path` CLI option.

## Emitting the checkpoints

During the simulation, one can enable checkpoint saving. The appropriate stanza in the input file is:

```yaml
# Parameters for the checkpoint generation.
checkpoints:
  # Whether it's enabled.
  enabled: true

  # Write new checkpoint every x time (in simulation time).
  save every: 1e2 tau

  # Format string for the checkpoint files. The singular parameter of the
  # format file is time (in the units of tau).
  path format: "ckpt/%.3f"
```