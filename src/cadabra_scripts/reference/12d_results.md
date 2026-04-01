Here we summarize the results of known results related to 12d interpretation of 10d type IIB effective actions

References: 1506.06757, 2507.07934

# Eight-derivative couplings from 12d

We follow [1506.06757] in deriving the 10d type IIB eight-derivative
couplings of gravity and axio-dilaton from a 12d perspective.

---

## Preliminaries

The two-derivative 10d effective action for the metric $g_{mn}$ and axio-dilaton
$\tau = C_{(0)} + i e^{-\Phi}$ is (in $\mathrm{SL}(2,\mathbb{R})$-invariant Einstein frame)

$$S^{(0)} = \int \bigl(R - 2\,P_m\,\bar{P}^m\bigr)\,{*_{10}}\,1\,$$

where $R$ is the 10d Ricci scalar and $P_m$ packages the $\tau$ kinetic term
in a way that makes the $\mathrm{SL}(2,\mathbb{R})/\mathrm{U}(1)$ coset structure manifest.

**Axio-dilaton quantities.**
Writing $\tau = \tau_1 + i\tau_2$, the complex vielbein, its conjugate,
and the composite U(1) connection are

$$P_m = \frac{i}{2\tau_2}\,\nabla_m\tau\,,\qquad \bar{P}_m = -\frac{i}{2\tau_2}\,\nabla_m\bar\tau\,,\qquad Q_m = -\frac{1}{2\tau_2}\,\nabla_m\tau_1 = \frac{i}{2}\bigl(P_m - \bar P_m\bigr). $$

Covariant derivatives (with respect to both the Levi-Civita and U(1) connections) are

$$D_m P_n = \nabla_m P_n - 2i\,Q_m\,P_n\,,\qquad D_m \bar P_n = \nabla_m \bar P_n + 2i\,Q_m\,\bar P_n\,. \tag{3}$$

The antisymmetrised derivative satisfies $D_{[m}P_{n]} = 0$ as long as $\tau$ is globally well-defined (no 7-brane sources).

---

## 12d metric embedding

The key F-theory idea is that the axio-dilaton $\tau$ may be identified
with the complex-structure modulus of an auxiliary torus.
Define a 12d metric $\hat{g}_{MN}$ with $M = 0,\dots,11$ by

$$\hat{g}_{MN} = \begin{pmatrix} g_{mn} & 0 & 0 \\[4pt] 0 & \dfrac{1}{\tau_2} & \dfrac{\tau_1}{\tau_2} \\[6pt] 0 & \dfrac{\tau_1}{\tau_2} & \dfrac{\tau_1^2+\tau_2^2}{\tau_2} \end{pmatrix}, \tag{4}$$

where $g_{mn}$ is the 10d Einstein-frame metric and the lower-right $2\times2$
block is the flat metric on a torus $T^2$ with modular parameter $\tau$.
The torus has fixed (unphysical) volume, and only 10d derivatives act:
$\partial_M = (\partial_m,\,0,\,0)$.

The two-derivative action (1) then lifts to

$$S^{(0)} = \int \hat{R}\,{*_{10}}\,1\,, \tag{5}$$

where $\hat R = \hat R_{MN}\hat g^{MN}$ is the 12d Ricci scalar.
There is no 12d promotion of the 10d measure.
From this one can show that $\hat{R} = R - 2 |\partial P|^2$

## 12d Riemann tensor

It had been given in [2507.07934] that according to the 12d metric ansatze, one can evalaute that in complexified coordinates,

$$
\hat{R}_{z\bar{z}z\bar{z}}
=-\frac{1}{4}|P|^2,\quad
\hat{R}_{mnz\bar{z}}=P_{[m}\bar{P}_{n]},\quad
R_{mzn\bar{z}} = -\frac{P_m\bar{P}_n}{2},\quad
R_{mznz}=\frac{D_m P_n}{2}
$$



## Non-holomorphic Eisenstein series

At the eight-derivative level, the known parity-even part of the 10d
effective action is proportional to $f_0\,t_8\,t_8\,R^4$.
The $\mathrm{SL}(2,\mathbb{Z})$-invariant, non-holomorphic Eisenstein series of weight $3/2$
that serves as the modular coefficient is

$$f_0(\tau,\bar\tau) = \sum_{(m,n)\neq(0,0)} \frac{\tau_2^{3/2}}{|m + n\tau|^3}\,. \tag{6}$$

For large $\tau_2$ (weak coupling, $g_s = e^{\langle\Phi\rangle} \ll 1$) this has the expansion

$$f_0(\tau,\bar\tau) = 2\zeta(3)\,\tau_2^{3/2} + \frac{2\pi^2}{3}\,\tau_2^{-1/2} + \mathcal{O}(e^{-\tau_2})\,, \tag{7}$$

where the first term is tree-level, the second is one-loop,
and the exponentials are D-instanton contributions.

---

## $\hat t_8\,\hat t_8\,\hat R^4$ in 12d

The $t_8$ tensor is defined as the product-of-metrics structure appearing
in the four-graviton amplitude (see [Minasian:2015bxa], appendix A):

$$\hat t_8^{N_1\dots N_8} = \tfrac{1}{16}\Bigl( -2\bigl(\hat g^{N_1N_3}\hat g^{N_2N_4}\hat g^{N_5N_7}\hat g^{N_6N_8} + \hat g^{N_1N_5}\hat g^{N_2N_6}\hat g^{N_3N_7}\hat g^{N_4N_8} + \hat g^{N_1N_7}\hat g^{N_2N_8}\hat g^{N_3N_5}\hat g^{N_4N_6}\bigr) \\ +8\bigl(\hat g^{N_2N_3}\hat g^{N_4N_5}\hat g^{N_6N_7}\hat g^{N_8N_1} + \hat g^{N_2N_5}\hat g^{N_6N_3}\hat g^{N_4N_7}\hat g^{N_8N_1} + \hat g^{N_2N_5}\hat g^{N_6N_7}\hat g^{N_8N_3}\hat g^{N_4N_1}\bigr) \\ -(N_1\leftrightarrow N_2) - (N_3\leftrightarrow N_4) - (N_5\leftrightarrow N_6) - (N_7\leftrightarrow N_8) \Bigr)\,. \tag{8}$$

The 12d promotion of $t_8\,t_8\,R^4$ is

$$S_{t_8}^{(3)} = \int f_0(\tau,\bar\tau)\;\hat t_8\,\hat t_8\,\hat R^4\;{*_{10}}\,1\,, \tag{9}$$

where

$$\hat t_8\,\hat t_8\,\hat R^4 = \hat t_8^{M_1\dots M_8}\,\hat t_{8\,N_1\dots N_8}\; \hat R^{N_1N_2}{}_{M_1M_2}\, \hat R^{N_3N_4}{}_{M_3M_4}\, \hat R^{N_5N_6}{}_{M_5M_6}\, \hat R^{N_7N_8}{}_{M_7M_8}\,. \tag{10}$$

**Goal.**
Expanding (10) using the 12d metric decomposition (4) yields the full 10d action
involving $R_{mnrs}$, $P_m$, $\bar P_m$, $D_m$, and $R_{mn}$.
This is equation (B.10) of [Minasian:2015bxa], which we aim to reproduce.
