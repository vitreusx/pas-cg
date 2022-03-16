# Internal units. Quantities

All the physical quantities in the input data can be given units. The quantity
format is either `[number]`, in which case it is said to be dimensionless,
or `[numerical value] [unit]`.

The list of predefined units is: "f77unit" (internal unit of length, equal to 5
Angstrem), "A" (angstrom), "nm", "m", "ns", "tau" (internal unit of time, equal
to 1 ns), "micros" (microseconds), "ms", "s", "atom", "mol", "eps" (internal
unit of energy, equal to $\approx 1.5 \textrm{kcal}/\textrm{mol}$), "kcal", "J"
, "kB" (Boltzmann constant), "K", "kg", "amu", "f77mass" (internal unit of mass,
equal to average residue mass), "e", "C", "Amp" (Ampere, since "A" is taken by
angstrom), "c", "H", "mu_0" (magnetic permittivity of the vacuum), "eps_0" (
electric permittivity of the vacuum), "rad", "deg".

The units can be divided or raised to nth power (with `**`) - in general, any
arithmetic expression involving these predefined units can be entered.

```{warning}
The system of units is inconsistent. From `f77unit`, `tau` and `eps` one could derive the base unit of mass, yet `f77mass` is **not** equal to it.
```

For more info,
see [`units.h`](https://github.com/vitreusx/pas-cg/blob/main/cg/include/cg/utils/units.h)
and [`quantity.cpp`](https://github.com/vitreusx/pas-cg/blob/main/cg/src/utils/quantity.cpp)
.