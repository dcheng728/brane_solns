---
layout: default
title: Implications of the dimensional reduction ansatz
---

# Implications of the dimensional reduction ansatz

## 11d to 10d ansatz

In 11d, we are measuring length in units of $l_P$, and in 10d we are using $l_S$, both of which has dimension of length.

We use 11d coordinates of $(\vec{x},u)$, where $x$ are Minkowski coordinates with dimension length and $u \sim u + 2\pi R$ is a coordinate with length dimension in the 11th dimension, whose length scale is $R$.

It is natural to write down a relation between the 11d and 10d metrics

$$
ds_{11}^2 = ds_{10}^2 + (du + A_m dx^m)^2
$$

with implicit unit conversions.
Because as numerical values $ds_{10}^2$ is measured as quantities of $l_S^2$, while $ds_{11}^2$ is measured in $l_P^2$, and $du$ is measured in $R$.
To make the unit conversions explicit, we divide dimensionful displacements such as $dx, ds^2, du$ by the natural length units they are measured in, to obtain a dimensionless equation, for numerical measurereadings that one would see on a ruler

$$
\frac{ds_{11}^2}{l_P^2}
=\frac{l_S^2}{l_P^2}
\left[
    \frac{ds_{(S),10}^2}{l_S^2}+
    \left(\frac{R}{l_S}\right)^2\left(\frac{du}{R} + \frac{\vec{A}}{R/l_S}\cdot \frac{d\vec{x}}{l_S}\right)^2
\right]
$$

We can obtain a dimensionless equation by dividing every dimensionful quan

The most general dimensional reduction formula, with the correct units is

$$
\frac{ds_{11}^2}{l_P^2}
=\frac{l_S^2}{l_P^2}\left[
    \frac{ds_{(S),10}^2}{l_S^2}+\left(\vec{A}\cdot \frac{d\vec{x}}{l_S}\right)^2\right]
    +\frac{R^2 d\theta^2}{l_P^2}
    +2\frac{l_S}{l_P}\frac{Rd\theta}{l_P}\left(\vec{A}\cdot \frac{d\vec{x}}{l_S}\right)
$$

In 10d string frame metric $ds_{(S),10}^2$ and $d\vec{x}$ are most naturally measured in $l_S$.
The compact circle $R d\theta$ originally came from 11d so it is naturally measured in $l_P$.
Then the necessary factors were added from dimensional analysis.

We can extract a dimensionless quantity $R/l_S$ and redefine $\frac{\vec{A}}{R/l_S} \to \vec{A}$, then the metric is

$$
\frac{ds_{11}^2}{l_P^2}
=\frac{l_S^2}{l_P^2}
\left[
    \frac{ds_{(S),10}^2}{l_S^2}
    +\frac{R^2}{l_S^2}\left(d\theta+\vec{A}\cdot \frac{d\vec{x}}{l_S}\right)^2
\right].
$$

Letting it be implicit that $ds_{11}^2$ is measured in $l_P^2$ and $ds_{10}^2$ is measured in $l_S^2, we have

$$
ds_{11}^2 = \frac{l_S^2}{l_P^2}\left[
    ds_{(S),10}^2 + \frac{R^2}{l_P^2}(d\theta+\vec{A}\cdot d\vec{x})^2
\right].
$$

Substituting in the standard KK reduction ansatz and replacing $e^\Phi$ with its vacuum expectation $g_s$, we find the relations

$$
l_P^3 = g_s l_S^3,\quad
R = g_s l_S.
$$

## Torus setup

## 12d to 10d ansatz

Similar to the above analysis, we write down a general form for 12d spacetime with length unit $l_F$ reduced diagonally on a torus to 10d with length unit $l_S$:

$$
\frac{ds_{12}^2}{l_F^2}
=\frac{l_S^2}{l_F^2}
\left[
    \frac{ds_{(S),10}^2}{l_S^2} + \frac{\sqrt{M_s}}{l_S^2\tau_2}
    \left[(du+\tau_1 dv)^2 + \tau_2^2 dv^2\right]
\right].
$$

We are adopting to the same notation in [12d_unification_incorrect.md](12d_unification_incorrect.md), where $u,v$ have dimension of length and the torus has volume $4\pi^2 \sqrt{M_s}l_S^2$. One can also choose to measure the torus with $l_F$, in which case the torus has volume $4\pi^2\sqrt{M_F}l_F^2$ , and the two choices are related by $\sqrt{M_s}l_S^2 = \sqrt{M_F}l_F^2$.


The 12d reduction ansatz to 10d written in string frame is

$$
ds_{12}^2 = e^{-\Phi/2}ds_{(S),10}^2 + e^\Phi\left[(du+Cdv)^2 + e^{-2\Phi}dv^2\right]
$$

The $ds_{12}^2$ must be measured in the 12d length unit $l_F$ while the string frame metric is measured in $l_S$. The coordinates $du, dv$ are measured in $l_S$ had we chosen to work with $\sqrt{M_s}$.
Putting the units back

$$
\frac{ds_{12}^2}{l_F^2} = e^{-\Phi/2}
\left[
    \frac{ds_{(S),10}^2}{l_S^2} + \frac{e^{\Phi/2}}{l_S^2e^{-\Phi}}\left[(du+Cdv)^2 + e^{-2\Phi}dv^2\right]
\right].
$$

Substituting $e^\Phi$ with $g_s$, we obtain the relations

$$
l_F^4 = g_s l_S^4,\quad
R = g_s^{1/4}l_S = l_F
$$

where $R^2 = \sqrt{M_S}l_S^2 = \sqrt{M_F}l_F^2$ gives the volume of the torus $4\pi^2 R^2$.
So the compactification is done on a torus whose volume coincides with the 12d length scale.