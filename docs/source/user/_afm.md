# Simulations with Atomic Force Microscope

## Parameter file entry

```yaml
# Options for pulling with Atomic Force Microscope tips.
AFM simulations:
  # Whether to perform such simulations.
  perform: false

  # Type of the AFM simulation. Options:
  # - "stretch": after equilibration, the first and last residues of the first
  #   chain are pulled apart. To be more specific, the first residue is kept in
  #   place by a not moving AFM tip of the "const velocity" type, and the last
  #   one is pulled in the opposite direction by an AFM tip of either "const
  #   velocity" or "const force" type.
  # - "pull, then release": before equilibration, the last residue of the first
  #   chain is pulled by an AFM tip in a direction from the first to the last
  #   residue.
  type: stretch

  # Time window, over which to average the force/velocity of the residues. See
  # the discussion above about the precise meaning of the "averaging".
  average stats over: 100 tau

  # Parameters for the "pull, then release" simulation type.
  pull-release params:
    # Time over which to pull the
    time: 100 tau

  # Parameters for the "AFM tips".
  tip params:
    # Type of the AFM tip used on the last residue of the chain. Options:
    # - "const velocity": the AFM tip is attached to the residue with a harmonic
    #   force, and moves in the prescribed direction.
    # - "const force": the "AFM tip" exerts a constant force on the residue,
    #   so in effect it's not an AFM tip as much as it's just an added force to
    #   the residue.
    type: const velocity

    # Parameters for the "const velocity" option.
    velocity params:
      # Velocity of the AFM tip.
      velocity: 5e-3 A/tau

      # H1 parameter in the harmonic force V = H1 d^2 + H2 d^4.
      H1: 30.0 eps/A**2

      # H2 parameter in the harmonic force V = H1 d^2 + H2 d^4.
      H2: 0.0 eps/A**2

    # Parameters for the "const force" option.
    force params:
      # Force of the "AFM tip".
      force: 1e-2 eps/(A*tau)
```