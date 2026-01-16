# Type IIB supergravity action

These are given in [^PolchinskiVol2]

## String frame (split into NS/R/CS)

$$
\begin{equation}
S_{\mathrm{IIB}}=S_{\mathrm{NS}}+S_{\mathrm{R}}+S_{\mathrm{CS}}.
\end{equation}
$$

$$
\begin{equation}
S_{\mathrm{NS}}=\frac{1}{2\kappa_{10}^2}\int d^{10}x(-G)^{1/2}e^{-2\Phi}\left(R+4\partial_{\mu}\Phi\partial^{\mu}\Phi-\frac{1}{2}|H_3|^2\right).
\end{equation}
$$

$$
\begin{equation}
S_{\mathrm{R}}=-\frac{1}{4\kappa_{10}^2}\int d^{10}x(-G)^{1/2}\left(|F_1|^2+|\tilde F_3|^2+\frac{1}{2}|\tilde F_5|^2\right).
\end{equation}
$$

$$
\begin{equation}
S_{\mathrm{CS}}=-\frac{1}{4\kappa_{10}^2}\int C_4\wedge H_3\wedge F_3.
\end{equation}
$$

## Einstein frame (directly from string frame)

We define the Einstein-frame metric by the Weyl rescaling

$$
\begin{equation}
G^{(E)}_{\mu\nu}=e^{-\Phi/2}G_{\mu\nu}.
\end{equation}
$$

In terms of $G_E$ the action becomes

$$
\begin{equation}
S_{\mathrm{IIB}}=S^{(E)}_{\mathrm{NS}}+S^{(E)}_{\mathrm{R}}+S_{\mathrm{CS}},
\end{equation}
$$

$$
\begin{equation}
S^{(E)}_{\mathrm{NS}}=\frac{1}{2\kappa_{10}^2}\int d^{10}x(-G_E)^{1/2}\left(R_E-\frac{1}{2}(\partial\Phi)^2-\frac{1}{2}e^{-\Phi}|H_3|^2\right),
\end{equation}
$$

$$
\begin{equation}
S^{(E)}_{\mathrm{R}}=-\frac{1}{4\kappa_{10}^2}\int d^{10}x(-G_E)^{1/2}\left(e^{2\Phi}|F_1|^2+e^{\Phi}|\tilde F_3|^2+\frac{1}{2}|\tilde F_5|^2\right).
\end{equation}
$$

## Einstein frame (SL(2)-covariant form)

$$
\begin{equation}
S_{\mathrm{IIB}}=\frac{1}{2\kappa_{10}^2}\int d^{10}x(-G_E)^{1/2}\left(R_E-\frac{\partial_{\mu}\bar\tau\partial^{\mu}\tau}{2(\mathrm{Im}\tau)^2}-\frac{\mathcal{M}_{ij}}{2}\tilde F_3^i\cdot\tilde F_3^j-\frac{1}{4}|\tilde F_5|^2\right)-\frac{\epsilon_{ij}}{8\kappa_{10}^2}\int C_4\wedge\tilde F_3^i\wedge\tilde F_3^j.
\end{equation}
$$

## References

[^PolchinskiVol2]: J. Polchinski, *String Theory, Volume 2: Superstring Theory and Beyond*, Cambridge University Press (1998).
