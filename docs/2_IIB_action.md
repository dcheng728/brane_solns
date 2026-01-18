---
layout: default
title: 2. Type IIB supergravity action
---

# Type IIB supergravity action

These are given in [^Polchinski98]

## String frame (split into NS/R/CS)

$$
\begin{equation}
S_{\mathrm{IIB}}=S_{\mathrm{NS}}+S_{\mathrm{R}}+S_{\mathrm{CS}}.
\end{equation}
$$

$$
\begin{equation}
S_{\mathrm{NS}}=\frac{1}{2\kappa_{10}^2}\int d^{10}x\sqrt{-g}e^{-2\Phi}\left(R+4\partial_{\mu}\Phi\partial^{\mu}\Phi-\frac{1}{2}|H_3|^2\right).
\end{equation}
$$

$$
\begin{equation}
S_{\mathrm{R}}=-\frac{1}{4\kappa_{10}^2}\int d^{10}x\sqrt{-g}\left(|F_1|^2+|\tilde F_3|^2+\frac{1}{2}|\tilde F_5|^2\right).
\end{equation}
$$

$$
\begin{equation}
S_{\mathrm{CS}}=-\frac{1}{4\kappa_{10}^2}\int C_4\wedge H_3\wedge F_3.
\end{equation}
$$

## Einstein frame 

The Einstein frame is the frame which the Einstein-Hilbert action takes the canonical form $\sqrt{-g} R$.
However this only determines the metric upto a constant factor

$$
\begin{equation}
g_{mn}^{(E)}
\equiv e^{-\frac{\Phi - \Phi_0}{2}}g^{(S)}_{mn}.
\end{equation}
$$

Then one finds that

$$
\begin{equation}
\sqrt{-g}\bigg|_{g_{mn}^{(S)}} = e^{\frac{5}{2}(\Phi-\Phi_0)}\sqrt{-g}\bigg|_{g_{mn}^{(E)}},\quad
(\partial\Phi)^2\bigg|_{g_{mn}^{(S)}} = e^{\frac{(\Phi_0-\Phi)}{2}}(\partial\Phi)^2\bigg|_{g_{mn}^{(E)}},\quad
|F_p|^2\bigg|_{g_{mn}^{(S)}} = e^{\frac{p(\Phi_0-\Phi)}{2}}|F_p|^2\bigg|_{g_{mn}^{(E)}},
\end{equation}
$$

as well as[^Bento23]

$$
\begin{equation}
\sqrt{-g} e^{-2\Phi}R\bigg|_{g^{(S)}_{mn}}
=\frac{1}{e^{2\Phi_0}}\sqrt{-g}\left[R - \frac{9}{2}(\partial\Phi)^2\right]\bigg|_{g^{(E)}_{mn}}.
\end{equation}
$$

Substituting these into the string frame action and letting $e^{\Phi_0}$ be denoted by $g_s$, we find the Einstein frame action

$$
\begin{equation}
\begin{aligned}
2\kappa_{10}^2g_s^2S_{IIB}
&=\int d^{10}x\sqrt{-g}\left[
R-\frac{1}{2}\frac{\partial\tau\partial\bar\tau}{\tau_2^2}
\right]\\
&-\frac{g_s}{2}\int d^{10}x\sqrt{-g}\left[
e^{-\Phi}|H_3|^2+e^\Phi|\tilde{F}_3|^2
\right]\\
&-g_s^2\left[\frac{1}{4}\int d^{10}x\sqrt{-g}
|\tilde F_5|^2
+\frac{1}{2}\int C_4\wedge H_3\wedge F_3\right].
\end{aligned}
\end{equation}
$$

We note that above we adopted the convention where $\Phi_0 \neq 0$, below is an alternative convention.

### Alternative convention (Polchinski)

An alternative convention can be found in Polchinski[^Polchinski98], which we present as

$$
\begin{equation}
g^{(E)}_{mn}=e^{-\Phi^{(E)}/2}g_{mn}^{(S)},\quad
\Phi^{(E)}\equiv \Phi - \Phi_0,
\end{equation}
$$

where $\Phi$ is the dilaton in string frame.

In terms of $g^{(E)}_{mn}$ the action becomes


$$
\begin{equation}
\begin{aligned}
S_{\mathrm{IIB}}&=\frac{1}{2\kappa_{10}^2}\int d^{10}x\sqrt{-g}
\Bigg[\left(R-\frac{\partial_{\mu}\bar\tau\partial^{\mu}\tau}{2\tau_2^2}\right)-\frac{1}{2}\left(e^{-\Phi^{(E)}}|H_3|^2+e^{\Phi^{(E)}}|\tilde F_3|^2\right)-\frac{1}{4}|\tilde F_5|^2\Bigg]\\
&\quad-\frac{1}{4\kappa_{10}^2}\int C_4\wedge H_3\wedge F_3.
\end{aligned}
\end{equation}
$$

### Comments

It appears that the different conventions of $g_s$ arises from the ambiguity in the definition of the Einstein frame metric, and that this ambiguity can be absorbed by redefining the dilaton in the Einstein frame.
So it appears that one has to decide between the two:

1. Accept explicit $g_s = e^{\Phi_0}$, then the non-compact scalar $\Phi$ of type IIB in string frame and Einstein frame are the same field
2. Set $g_s = 1$, and accept the string frame and Einstein frame non-compact scalar $\Phi$ to differ by a constant: $\Phi^{(E)} = \Phi - \Phi_0$.

<!-- ## Einstein frame (SL(2)-covariant form)

$$
\begin{equation}
S_{\mathrm{IIB}}=\frac{1}{2\kappa_{10}^2}\int d^{10}x(-G_E)^{1/2}\left(R_E-\frac{\partial_{\mu}\bar\tau\partial^{\mu}\tau}{2(\mathrm{Im}\tau)^2}-\frac{\mathcal{M}_{ij}}{2}\tilde F_3^i\cdot\tilde F_3^j-\frac{1}{4}|\tilde F_5|^2\right)-\frac{\epsilon_{ij}}{8\kappa_{10}^2}\int C_4\wedge\tilde F_3^i\wedge\tilde F_3^j.
\end{equation}
$$ -->

## References

[^Polchinski98]: J. Polchinski, *String Theory, Volume 2: Superstring Theory and Beyond*, Cambridge University Press (1998).

[^Bento23]: B. V. Bento, D. Chakraborty, S. Parameswaran, and I. Zavala, *A guide to frames, 2 $\pi$’s, scales and corrections in string compactifications* (2023), [arXiv:2301.0517 [hep-th]](https://arxiv.org/abs/2301.05178).
