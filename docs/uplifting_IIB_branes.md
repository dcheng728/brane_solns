# Uplifting IIB branes to 12d

## Type IIB branes

We have review the type IIB brane solutions supported by form-fields in [review](review.pdf). They are given by

with a single harmonic function $H$ on the transverse space, where for codimension $\ge 3$ one may take $H=1+\frac{Q}{r^{7-p}}$. In string frame, the elementary solutions are:

- **Fundamental string (F1)** (supported by $B_2$). With worldvolume coordinates $x^{0,1}$ and transverse coordinates $x^m$,

    $$
    \begin{equation}
    ds_{(S),10}^2=H^{-1}(-dt^2+dx_1^2)+dx^mdx^m,\qquad e^{2\Phi}=g_s^2H^{-1},\qquad B_{01}=H^{-1}-1.
    \end{equation}
    $$

- **NS5-brane** (supported by $H_3=dB_2$). With worldvolume coordinates $x^{0,\dots,5}$,

    $$
    \begin{equation}
    ds_{(S),10}^2=dx^\mu dx^\mu+H\,dx^mdx^m,\qquad e^{2\Phi}=g_s^2H,\qquad H_3=\star_4 dH.
    \end{equation}
    $$

- **D$p$-branes** (supported by RR fields), where for type IIB we have $p=-1,1,3,5,7,9$. With worldvolume coordinates $x^{0,\dots,p}$,

    $$
    \begin{equation}
    ds_{(S),10}^2=H^{-\frac{1}{2}}dx^\mu dx^\mu+H^{\frac{1}{2}}dx^mdx^m,\qquad e^{\Phi}=g_s\,H^{\frac{3-p}{4}},\qquad C_{0\cdots p}=H^{-1}-1.
    \end{equation}
    $$

## 12d to 10d reduction ansatz, assuming SO(12) -> SO(10)xSO(2)
What should be the reduction ansatz from 12 to 10?
If we assume that in 12d one should not know about the torus, hence know nothing about 10d, and that only once the torus is introduced the 10d theory is obtained via 12 = 10 + 2.
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
ds_{12}^2 = e^{-\frac{1}{2}\Phi}ds_{(S),10}^2 + e^\Phi[(du+Cdv)^2 + e^{-2\Phi}dv^2]
\end{equation}
$$

## 12D uplift of the various branes

Setting $C = 0$ here.

### D(-1)

### F1

### D1

### D3

### NS5

### D5

### D7