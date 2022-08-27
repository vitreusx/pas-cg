# Simulation box settings 

(squeezing)=
## Squeezing

(finding-force-minimum)=
## Finding force minimum

(oscillations)=
## Oscillations

(pulling-at-the-end)=
## Pulling at the end

## Parameter file stanza

```yaml
# Simulation box settings, as well as the parameters detailing simulation-box-
# related simulation protocols, such as squeezing, oscillations etc.
simulation box:
  # Initial size specification.
  initial size:
    # Type of the derivation of the initial simulation box size. Options:
    # - "sufficient": the box is such as to contain all the residues. See
    #   `sufficient box params` for more details.
    # - "from CRYST1 record": the box is set to one spanning from [0, 0, 0] to
    #   [a, b, c] where the latter values are taken from CRYST1 record in the
    #   PDB file. Note: this usually only makes sense if the periodic boundary
    #   conditions are enabled, or there are no walls.
    # - "keep from SAW": the box is set to one which was used in the SAW
    #   procedure with the periodic boundary conditions (Note: this is **not**
    #   the same as the CRYST1 box, as the latter is not actually used in the
    #   SAW procedure, even if periodic boundary conditions are used)
    # - "infinite": the simulation box is boundless; of course, this setting
    #   only makes sense with no walls.
    type: sufficient

    # Further parameters for the "sufficient" option for the "type" property
    # above.
    sufficient box params:
      # A multiple of the bond size (derived as mentioned previously), which
      # gets added in each direction to the minimal box containing the residues.
      pad (x bond): 2.0

      # Maximum density of the box. Optional. To be more specific, if set, the
      # simulation box necessarily contains a cube with a size derived so as to
      # make it have the specified molar density, in a fashion similar to an
      # option "start box" in the section about morphing into SAW, centered in
      # zero.
      maximum density: 1e-4 residue/A**3

      # Whether the generated simulation box should be cubic or not.
      cubic: false

  # Parameters for the simulation box walls.
  walls:
    # Types of the walls. This can be set separately for x/y/z axes, as well as
    # set for all axes at the same time.
    # The key options are:
    # - "x axis", "y axis", "z axis": separately;
    # - "all axes": all at once.
    # The wall types are:
    # - "void": no walls on the given axis/axes;
    # - "periodic": periodic boundary conditions on a given
    #   axis/axes (this is, one must admit, a mild misuse of the language, but
    #   putting it here seems like a most fitting option.)
    # - "solid": the walls being "solid", i.e. described by potential V ~ d^{-9}.
    # Note that, if enabled, the walls can change the type to an attractive one
    # during the simulation; this usually only makes sense if the type of the
    # wall is "solid".
    all axes: void

    # Threshold for the wall interactions. (Note: equivalent to `walmindist`)
    threshold: 5.0 A

    # Parameters for the "solid" wall type.
    solid wall params:
      # Depth of the potential.
      depth: 4.0 eps

    # Parameters for the LJ attractive wall type.
    lj wall params:
      # Depth of the potential.
      depth: 4.0 eps

      # Timespan, over which the target depth is acquired quasi-adiabatically.
      cycle duration: 0.0 tau

      # A multiplier for the minimum of the radius, above which the connection
      # of a residue with the wall is considered to be broken.
      breaking dist factor: 1.5

      # Connection limit for each wall.
      # Options:
      # - a number: a fixed number of connections.
      # - "auto": a value derived via a fairly complex formula which I won't
      #   reproduce here.
      connection limit: auto

    # Parameters for the harmonic attractive wall type.
    harmonic wall params:
      # The "depth": V = HH1 d^2
      HH1: 30.0 eps/A**2

      # Connection limit for each wall. See the description in the
      # "lj wlal params" for details.
      connection limit: auto

  # Common parameters for squeezing, oscillations etc.
  common:
    # A "resting period" between various phases of the simulation.
    rest period: 10e3 tau

    # Maximum (regular) velocity of wall movement, reached with constant
    # acceleration over a given distance.
    target velocity: 5e-3 A/tau

    # Distance, over which the walls reach the full/target velocity.
    distance before acquiring full velocity: 12.0 A

    # A span of time, over which to average the wall forces. The value of this
    # force average is not computed with relation to the current time point -
    # rather, simulation is divided into "time windows", and the current value of
    # the force average is equal to the average of forces in the *previous* window;
    # in the first window, the value is undefined.
    average forces over: 100 tau

  # Squeezing phase settings.
  squeezing:
    # Whether to perform box squeezing.
    perform: false

    # Density of the target box.
    target density: 1e-3 residue/A**3

    # Velocity used when the simulation box volume is twice above the target
    # volume - this being done to speed up the simulation, as one presumes that
    # one needs not be that careful when there are barely any residues
    # interacting with the walls.
    velocity above 2V: 2e-2 A/tau

  # "Finding force minimum" phase.
  finding force minimum:
    # Whether to perform the phase.
    perform: false

    # Force for max velocity. The velocity of the wall is proportional to the
    # forces acting on the wall, clipped at a given force threshold.
    force for max velocity: 0.4 eps/A

  # Oscillations phase.
  oscillations:
    # Whether to perform the phase.
    perform: false

    # Type of wall movement. Options:
    # - "squeeze": in which the z-axis walls are moved closer/apart.
    # - "shear": in which the z-axis walls are moved along the x-axis.
    type: squeeze

    # Number of oscillation cycles to perform.
    num of cycles: 6

    # Parameters for the amplitude of the oscillations.
    amplitude:
      # Variant of the derivation of the amplitude. Options:
      # - "absolute": the amplitude value is as given in "absolute value" entry
      #   below;
      # - "relative": the amplitude value is the distance between the z-axis
      #   walls multiplied by the value of "relative value" entry below.
      variant: absolute
      absolute value: 10.0 A
      relative value: 0.1

    # Angular frequency of the oscillations.
    angular frequency: 1e-4 rad/tau

  # "Pulling at the end" phase. [Not implemented]
  pulling at the end:
    # Whether to perform the phase.
    perform: false

  # Parameters for attractive walls. [Not implemented]
  attractive walls:
    # When to turn on the attractive walls. Options: ??
    when: never
    # Type of attractive wall: Options: "lj", "harmonic".
    type: lj
```