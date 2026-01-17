---
layout: default
title: Torus setup
---

# Torus setup

Reference ^Polchinski98

## Torus metric

A torus is $S^1 \times S^1$ with two coordinates $(\sigma^1,\sigma^2)\sim (\sigma^1,\sigma^2) + 2\pi R (m,n)$.
Using (Diff $\times$ Weyl), it is always possible to bring the metric to be flat $ds^2 = (d\sigma^{1})^2+(d\sigma^{2})^2$, but, depending on the torus to begin with, the periodicity condition in general will not remain as that of a rectangle, but become that of a rotated parallelgram

$$
\begin{equation}
(\sigma^1,\sigma^2)\sim (\sigma^1,\sigma^2)+2\pi m(u^1,u^2)+2\pi n(v^1,v^2)
\end{equation}
$$

Using further Diff $\times$ Weyl that keeps the metric in flat form, the rotation can be removed, ie $u^a = (1,0)$, and the periodicity condition is

$$
\begin{equation}
(\sigma^1,\sigma^2) \sim (\sigma^1,\sigma^2) + 2\pi (m+nv^1,n v^2)
\end{equation}
$$

We have been working with $\mathbb R^2$ and using orthonormal basis vectors $\vec{e}_1, \vec{e}_2$, but we can also go to $\mathbb C$, and use basis vectors $1, i$, then the equivalence relation is 

$$
\begin{equation}
\sigma^1 + i\sigma^2 \sim \sigma^1 + i\sigma^2 +  2\pi m + 2\pi n \tau, \quad \tau = v^1 + iv ^2
\end{equation}
$$

In $\mathbb R^2$, recall we had already used DiffxWeyl to set the metric as flat $ds^2 = (d\sigma^1)^2 + (d\sigma^2)^2$, in $\mathbb C$, this is $ds^2 = |d\sigma^1 + id\sigma^2|^2$ . 

In $\mathbb C$, we can not only use $1,i$ as basis vectors but also use $1,\tau$, as bases, as long as they are linearly independent, ie $Im \tau \neq 0 $, which we will demand. Using 1,$\tau$ as bases, the components $\sigma^1,\sigma^2$ will change, ie $\sigma^1 + i\sigma^2 = \tilde{\sigma}^1 + \tau \tilde \sigma^2$, and the benefit here is that we get back the rectangular periodicity condition

$$
\begin{equation}
\tilde \sigma^1 + \tau \tilde \sigma^2 \sim 
\tilde \sigma^1 + \tau \tilde \sigma^2 + 2\pi (m,\tau n)
\end{equation}
$$

Thus we have, for these new coordinates

$$
\begin{equation}
ds^2 = |d\sigma^1 + d\sigma^2|^2
=|d\tilde \sigma^1 + \tau d\tilde \sigma^2|^2
\end{equation}
$$

Now rename, and get 

$$
\begin{equation}
\boxed{
ds^2 = |d\sigma^1 + \tau d\sigma^2|^2,\quad
(\sigma^1,\sigma^2)\sim(\sigma^1 + 2\pi m ,\sigma^2+2\pi n)
}.
\end{equation}
$$

Here the periodicity condition of $\sigma \in [0,2\pi]$ implies that they are dimensionless angular coordinates.
As in, they are given as some others dimensionful coordinates divided by some length $R/\sqrt{\tau_2}$ , which is the dimensionful length scale of the torus.
To put the proper dimensions back, we write

$$
\begin{equation}
\boxed{
ds^2 = \frac{R^2}{\tau_2}|d\sigma^1 + \tau d\sigma^2|^2=|dy^1+\tau dy^1|^2, \quad
\vec{y} \equiv \frac{R}{\sqrt{\tau_2}} \vec{\sigma} \in \left[0,\frac{2\pi R}{\sqrt{\tau_2}}\right],\quad
(y^1,y^2)
\sim \left(y^1 +  \frac{2\pi mR}{\sqrt{\tau_2}},y^2+ \frac{2\pi n R}{\sqrt{\tau_2}}\right)
}.
\end{equation}
$$

The volume of the torus is the determinant of the metric which is

$$
\begin{equation}
Vol(T_2) = 4\pi^2 R^2
\end{equation}
$$


## SL(2, Z)

In the Quantum Theory, path-integrating over all metrics on the torus becomes path-integrating over all in-equivalent values of $\tau$. As a general number in $\mathbb C$, different values of $\tau$ can be parameterizing the same torus, we would like to identify to what degree this is happening, we do not want to integrate over redundancy in our path integral, as we would only like to integrate over indeed physically different torus. To do so we look at the equivalence relation given, and identify what transformation on the torus parameter $\tau$ lets us stay within this equivalence class, ie take us from $\sigma^a$ to $\sigma^{\prime a}$ such that $[\sigma^a] = [\sigma^{\prime a}]$, there are two types of those transformations, correspond to two operations on the lattice.

- $\tau \to \tau + 1$ is equivalent to $\tau \to \tau $, $(m,n) \to (m+n,n)$
- $\tau \to 1/\tau$ is equivalent to $\tau \to \tau$, $(m,n) \to (n,m)$, followed by a constant scaling of the coordinates $\sigma^a \to \tau \sigma^a$ which leaves the metric invariant. 


These two operations generate the SL(2,Z) transformations on $\tau$: $\tau \to \frac{a\tau +b}{c^\tau + d}$, $ad-bc = 1$. So this is how one should think of the SL(2,Z) symmetry of the torus: it is a footprint of the periodic identification one has to make, when putting coordinates on $S^1 \times S^1$. In other words, whenever one writes the metric inside the torus as $ds^2 = |d\sigma^1 + \tau d\sigma^2|^2$, there must be this SL(2,Z) symmetry on $\tau$ because before the coordinate transformations, when the metric was written in a less agreeable form, there was a periodic ientification on $\sigma^1,\sigma^2$.
In the IIB literature, one often packages $\tau = C + ie^{-\phi}$, then the metric takes the form

$$
    e^\phi
    \begin{pmatrix}
        e^{-2\phi}+C^2 & C \\ C & 1
    \end{pmatrix}
$$

## References

[^Polchinski98]: J. Polchinski, *String Theory, Volume 1: An introduction to the bosonic string*, Cambridge University Press (1998).