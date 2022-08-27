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

## Parameter file entry

```yaml
# Parameters for the Langevin integrator.
langevin:
  # Whether it's enabled. Note: if no integrator is enabled, of which there is
  # only one at the moment, the time doesn't move forward.
  enabled: true

  # Gamma factor in the equation mr'' = -dV/dr - gamma r' + noise, where the
  # noise has variance sigma^2 = 2 gamma k_B T
  gamma factor: 2.0f 1/tau

  # Temperature; see the variance of the noise term above.
  temperature: 0.35f eps/kB

  # Integrator time step.
  dt: 5e-3f tau
```