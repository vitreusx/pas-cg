# Native dihedral angle FFs

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

## Parameter file format

```yaml
complex native dihedrals:
  enabled: boolean
  CDA: quantity [E]
  CDB: quantity [E]

simple native dihedrals:
  enabled: boolean
  CDH: quantity [E/Angle**2]
```

Only one of these two can be enabled at the same time.