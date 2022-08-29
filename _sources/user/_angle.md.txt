# Angle potentials

## Native bond angle FF

The native bond angle force field is described by:

$$V = \sum \frac{1}{2}k(\theta - \theta_0)^2$$

where the summation is over the triples with the native value of the bond angle
present.

The bond angle itself is determined as follows:

$$\theta_{123} = \arccos
\left[-\frac{\langle \vec{r_{12}}, \vec{r_{23}} \rangle}{r_{12} r_{23}}\right]$$

## Heurestic bond angle FF

The heurestic bond angle force field is described by:

$$ V = \sum f(\theta_i, \mathtt{type}_i)\\ f(\theta, \mathtt{type}) = \sum_
{i=0}^D c_i[\mathtt{type}] \theta^D $$

where the summation is over the triples with the native value of the bond angle
not present. In other words, we approximate the value of the potential by a set
of polynomials, one for each different _type_ of the triple of residues. The
type of the triple is derived in the following fashion:

1. First, we class the residues into glycines (code `G`), prolines (code `P`)
   and all other amino acids into a single category (code `X`);
2. Then, the type of the triple of residues corresponds to the pair consisting
   of the class of the second and of the third residue. Thus for example, if we
   have a triple of residues which are glycine, proline and valine, the type of
   the residue is `PX` (`P` for the proline, `X` for the valine).
3. Thus, there are nine total types: `GG`, `GP`, `GX`, ..., `XP`, `XX`.

## Native dihedral angle FFs

There are two variants of the native dihedral angle force field: a _simple_
variant and a _complex_ variant. Both of them are described by:

$$ V = \sum f(\varphi - \varphi_0)
$$

where the summation is over the quadruples with the native value of the dihedral
angle present. They differ in the choice of the function $f$. For the complex
variant, the choice is:

$$ f_{\textrm{comp}}(\Delta\varphi) = \mathtt{CDA}(1-\cos{\Delta\varphi}) +
\mathtt{CDB}(1-\cos{3\Delta\varphi})
$$

whereas for the simple variant, it is:

$$ f_{\textrm{simp}}(\Delta\varphi) = \frac{1}{2}\mathtt{CDH} {\Delta\varphi}^2
$$

The value of the dihedral angle is computed in a following fashion:

$$ \varphi_{1234}' =
\arccos\left[\frac{\langle \vec{r_{12}} \times \vec{r_{23}}, \vec{r_{23}} \times \vec{r_{34}} \rangle}{|\vec{r_{12}} \times \vec{r_{23}}| |\vec{r_{23}} \times \vec{r_{34}}|} \right]
\\ \varphi_{1234} = \begin{cases} \varphi_{1234}' & \langle \vec{r_{12}} \times
\vec{r_{23}}, \vec{r_{34}} \rangle \geq 0\\ -\varphi_{1234}' &
\mathrm{otherwise} \end{cases} $$

## Heurestic dihedral angle FF

The heurestic dihedral angle force field is described by:

$$ V = \sum f(\varphi_i, \mathtt{type}_i)\\ f(\varphi, \mathtt{type}) = a + b
\sin\varphi + c \cos\varphi + d\sin^2\varphi + e\cos^2\varphi + f\sin\varphi
\cos\varphi $$

where the summation is over the quadruples without the native value and the
coefficients are different for each $\mathtt{type}$. The type of the quadruple
is determined in the same fashion as with the heurestic bond angle types, with the
second and the third of the four residues determining the type.

## Parameter file stanza

```yaml
# Parameters for the (bond and dihedral) angle potentials of various types.
angle potentials:
  # Whether all the potentials (native bond angle, heurestic bond angle,
  # native dihedral angle, heurestic dihedral angle) are enabled. This can be
  # overriden in the force-specific stanzas below.
  all enabled: true

  # Whether the potentials for the dihedral angles are enabled.
  dihedral potentials enabled: true

  # Parameters for the heurestic bond angle potential.
  heurestic bond angles params:
    # Coefficients of the polynomial describing the potential, depending on the
    # type of residue (GLY, PRO and rest as X).
    coefficients: |
      type1,type2,a0,a1,a2,a3,a4,a5,a6
      G,G,20872.75597,-63260.52963,79318.301,-52680.17088,19554.21897,-3847.670279,313.6716916
      G,P,8222.83155,-25178.89003,31841.70634,-21290.04519,7941.501624,-1567.87285,128.0761621
      G,X,20872.75597,-63260.52963,79318.301,-52680.17088,19554.21897,-3847.670279,313.6716916
      P,G,34646.70029,-109957.2324,144423.9036,-100525.9874,39127.31928,-8079.214542,691.8417699
      P,P,10744.12043,-34148.94233,44818.66284,-31110.67875,12060.69185,-2479.723349,211.6367439
      P,X,34646.70029,-109957.2324,144423.9036,-100525.9874,39127.31928,-8079.214542,691.8417699
      X,G,15883.02041,-48923.1471,62360.6974,-42110.86572,15891.78309,-3178.490602,263.2916319
      X,P,16912.27207,-53570.09757,70150.19389,-48602.41198,18791.04978,-3844.690523,325.3085829
      X,X,15883.02041,-48923.1471,62360.6974,-42110.86572,15891.78309,-3178.490602,263.2916319

  # Parameters for the native bond angle potential.
  native bond angles params:
    # `CBA` value in the potential V = CBA (theta - theta_n)^2.
    CBA: 30.0 eps/rad**2

  # Parameters for the heurestic dihedral angle potential.
  heurestic dihedral angles params:
    # Coefficients for the equation describing the potential, depending on the
    # type of residue (see above).
    coefficients: |
      type2,type3,const,sin,cos,sin2,cos2,sin_cos
      G,G,0.133672207,-0.007576316,0.003786907,-0.124627339,0.425373566,-0.060606303,
      G,P,0.935285048,0.928786712,-0.18516837,0.015857805,0.2861951,0.072728001
      G,X,0.210489196,-0.00606094,0.202709724,-0.160512736,0.461339767,0.13333598
      P,G,0.233402207,-0.101516187,0.109235732,0.14906496,0.151803103,-0.742423775
      P,P,1.810497634,1.171212546,0.091084321,-0.254152251,0.557619284,-1.569694253
      P,X,0.661379307,0.115151794,0.428904959,0.200723546,0.100490651,-0.803028162
      X,G,0.198889776,0.018181703,-0.070746181,0.122076238,0.178719533,-0.624241103
      X,P,1.254229713,0.739393723,0.686217752,0.219188398,0.083115678,-0.790909154
      X,X,0.275933756,0.00606063,0.257226522,0.15460117,0.146208844,-0.448484074

  # Parameters for the native dihedral angle potential.
  native dihedral angles params:
    # Variant of the potential. Options: "complex", "simple".
    variant: complex

    # Complex variant parameters.
    complex variant params:
      # CDA parameter in the formula which I won't reproduce here.
      CDA: 0.66 eps

      # CDB parameter in the formula which I won't reproduce here.
      CDB: 0.66 eps

    # Simple variant parameters.
    simple variant params:
      # CDH parameter in the formula which I won't reproduce here.
      CDH: 3.33 eps/rad**2
```