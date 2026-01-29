---
layout: default
title: 3. Uplifting IIB branes to 12d
---

# Uplifting IIB branes to 12d

## Type IIB branes

We use the convention where $g^{(S)}_{mn} = e^{(\Phi-\langle\Phi\rangle)/2}g^{(E)}_{mn}$

We have reviewed the standard construction of type IIB brane solutions supported by form-fields in [review](1_review.md). They are given by a single harmonic function $H$ on the transverse space, where one may take $H=1+\frac{Q}{r^{7-p}}$. In string frame, the elementary solutions are:

- Fundamental string (F1): (supported by $B_2$). With worldvolume coordinates $x^{0,1}$ and transverse coordinates $x^m$,

$$
ds_{(S),10}^2=H^{-1}(-dt^2+dx_1^2)+dx^mdx^m,\qquad \frac{e^{\Phi}}{g_s}=H^{-1/2},\qquad
(B_2)_{01} = H^{-1}-1.
$$

- NS5-brane (supported by $H_3=dB_2$). With worldvolume coordinates $x^{0,\dots,5}$,

$$
ds_{(S),10}^2=dx^\mu dx^\mu+H\cdot dx^mdx^m,\qquad \frac{e^{\Phi}}{g_s}=H^{1/2},\qquad H_3=\star_4 dH.
$$

- D $p$ -branes (supported by RR fields), where for type IIB we have $p=-1,1,3,5,7,9$. With worldvolume coordinates $x^{0,\dots,p}$,

$$
ds_{(S),10}^2=H^{-\frac{1}{2}}dx^\mu dx^\mu+H^{\frac{1}{2}}dx^mdx^m,\qquad \frac{e^\Phi}{g_s}=H^{\frac{3-p}{4}},\qquad C_{0\cdots p}=H^{-1}-1.
$$

## 12d to 10d reduction ansatz, assuming SO(12) -> SO(10)xSO(2)
What should be the reduction ansatz from 12 to 10?
If we assume that in 12d one should not know about the torus, hence know nothing about 10d, and that only once the torus is introduced the 10d theory is obtained via the symmetry breaking 12 = 10 + 2.
The reduction ansatz should only break symmetry as 12 = 10 + 2.

$$
\begin{equation}
ds_{12}^2 = e^{f\Phi} ds_{(S),10}^2 + e^{g\Phi} ds_2^2
\end{equation}
$$

and $e^\Phi$ is the 10d dilaton for the shape moduli of the torus. 
There is no other non-compact scalar to use for comformal factors.
The metric $ds_2^2$ is the metric on the torus of constant volume, given by $e^\Phi, C$.

From F-theory, and the pp-wave uplift of the D(-1), we are led to $f = -\frac{1}{2}, g = 0$. In other words, we can write

$$
\begin{equation}
ds_{12}^2 = g_s^{1/2}e^{-\Phi/2}ds_{(S),10}^2 + e^\Phi[(dy+Cdt)^2 + e^{-2\Phi}dt^2]
\end{equation}
$$

## 12D uplift of the various branes

Setting $C = 0$ here.

### D(-1)

PP-wave

### F1

$$
\begin{equation}
ds_{12}^2=
H^{-3/4}ds_{1,1}^2 + H^{1/4}ds_8^2+g_sH^{-1/2}dz_1^2 + g_s^{-1}H^{1/2}dz_2^2,
\end{equation}
$$

$$
(\mathcal C_3)_{012} = H^{-1}-1.
$$

### D1

$$
\begin{equation}
ds_{12}^2 = H^{-3/4}ds_{1,1}^2 + H^{1/4}ds_8^2 + g_sH^{1/2}dz_1^2 + g_s^{-1}H^{-1/2}dz_2^2,
\end{equation}
$$

$$
(\mathcal C_3)_{013} = H^{-1}-1.
$$

### D3

$$
ds_{12}^2=
H^{-\frac{1}{2}}ds_{1,3}^2
+H^{\frac{1}{2}}ds_6^2
+g_s\cdot du^2+g_s^{-1}\cdot dv^2.
$$

### NS5

$$
ds_{12}^2 = 
H^{-1/4}ds_{1,5}^2 + H^{3/4}ds_4^2 + g_sH^{1/2}dz_1^2 + g_s^{-1}H^{-1/2}dz_2^2,
$$

$$
(\mathcal F_4)_{m_1 m_2 m_3 2} = \varepsilon_{m_1 m_2 m_3}{}^n\partial_n H
$$


### D5

$$
ds_{12}^2 =
H^{-1/4}ds_{1,5}^2 + H^{3/4}ds_4^2 + g_sH^{-1/2}dz_1^2 + g_s^{-1}H^{1/2}dz_2^2
$$

$$
(\mathcal F_4)_{m_1 m_2 m_3 3} = \varepsilon_{m_1 m_2 m_3}{}^n\partial_n H
$$

### D7

KK-monopole