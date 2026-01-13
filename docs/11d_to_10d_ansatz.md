# 11d to 10d ansatz

In 11d, we are measuring length with $l_p$ and in 10d we are using $l_s$, both of which has dimension of length.

Restoring units in $ds^2$, we use 11d coordinates of $(\vec{x},\theta)$, where $x$ are Minkowski coordinates with dimension length and $u$ is a dimensionless compact coordinate on a circle of radius $R$.
The most general dimensional reduction formula, with the correct units is

$$
\frac{ds_{11}^2}{l_p^2}
=\frac{l_s^2}{l_p^2}\left[
    \frac{ds_{(S),10}^2}{l_s^2}+\left(\vec{A}\cdot \frac{d\vec{x}}{l_s}\right)^2\right]
    +\frac{R^2 d\theta^2}{l_p^2}
    +2\frac{l_s}{l_p}\frac{Rd\theta}{l_p}\left(\vec{A}\cdot \frac{d\vec{x}}{l_s}\right)
$$

In 10d string frame metric $ds_{(S),10}^2$ and $d\vec{x}$ must be measured in $l_s$ since $l_s$ is the length unit in 10d.
The compact circle originally came from 11d coordinates so it is originally measured in $l_p$.

We can extract a dimensionless quantity $R/l_s$ and redefine $\frac{\vec{A}}{R/l_s} \to \vec{A}$, then the metric is

$$
\frac{ds_{11}^2}{l_p^2}
=\frac{l_s^2}{l_p^2}
\left[
    \frac{ds_{(S),10}^2}{l_s^2}
    +\frac{R^2}{l_s^2}\left(du+\vec{A}\cdot \frac{d\vec{x}}{l_s}\right)^2
\right].
$$

Letting it be implicit that $ds_{11}^2$ is measured in $l_p^2$ and $ds_{10}^2$ is measured in $l_s^2, we have

$$
ds_{11}^2 = \frac{l_s^2}{l_p^2}\left[
    ds_{(S),10}^2 + \frac{R^2}{l_p^2}(du+\vec{A}\cdot d\vec{x})^2
\right].
$$

Substituting in $l_p^3 = e^\Phi l_s^3$ and $R = e^\Phi l_s$ we have the standard KK reduction ansatz.