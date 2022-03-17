# Langevin PC integrator

The equations of motions are:

$$ m_i \vec{a} = \vec{F} - \gamma \vec{v} + \vec{\Gamma} $$

where $\gamma$ is the damping coefficient, and $\vec{\Gamma}$ is the thermal
white noise with variance $\sigma^2 = 2\gamma k_BT$. For the purposes of
integration, six extra vectors $\vec{y}^{(0)}, \ldots, \vec{y}^{(5)}$ are stored
per residue, which represent consecutive terms of the Taylor expansion:

$$ \vec{y}^{(k)} = \frac{{\Delta t}^k}{k!}\frac{d^k\vec{r}}{dt^k} $$

The update step consists of two parts:

1. The predictor part.

$$ \begin{align*} \varepsilon &\leftarrow \vec{y}^{(2)} - \frac{{\Delta t}^2}{2!
} \vec{a}\\ \vec{y}^{(0)} &\leftarrow \vec{y}^{(0)} - 3\varepsilon/16\\
\vec{y}^{(1)}&\leftarrow \vec{y}^{(1)} - 251\varepsilon/360\\ \vec{y}^{(2)}
&\leftarrow \vec{y}^{(2)} - \varepsilon\\ \vec{y}^{(3)} &\leftarrow \vec{y}^{(3)
}- 11\varepsilon/18\\ \vec{y}^{(4)} &\leftarrow \vec{y}^{(4)} - \varepsilon/6\\
\vec{y}^{(5)} &\leftarrow \vec{y}^{(5)} - \varepsilon/60 \end{align*} $$

2. The corrector part.

$$ \begin{align*} \vec{y}^{(0)} &\leftarrow \vec{y}^{(0)} + \vec{y}^{(1)} +
\vec{y}^{(2)} + \vec{y}^{(3)} + \vec{y}^{(4)} + \vec{y}^{(5)}\\ \vec{y}^{(1)}
&\leftarrow \vec{y}^{(1)} + 2\vec{y}^{(2)} + 3\vec{y}^{(3)} + 4\vec{y}^{(4)} +
5\vec{y}^{(5)
}\\ \vec{y}^{(2)} &\leftarrow \vec{y}^{(2)} + 3\vec{y}^{(3)} + 6\vec{y}^{(4)} +
10\vec{y}^{(5)}\\ \vec{y}^{(3)} &\leftarrow \vec{y}^{(3)} + 4\vec{y}^{(4)} +
10\vec{y}^{(5)}\\ \vec{y}^{(4)} &\leftarrow \vec{y}^{(4)} + 5\vec{y}^{(5)}
\end{align*}$$

```{warning}
In the Fortran program, the damping was added directly to $\vec{v}$; it's unclear whether this was intended or not.
```

## Parameter file entry

```yaml
langevin:
  enabled: boolean
  type: one of [normal, legacy]
  gamma factor: quantity [1/T]
  temperature: a quantity or a pair of quantities, in [Temp]
  dt: quantity [T]
```

If a pair of temperature quantities is provided, then the temperature is
varied (linearly) for different trajectories.