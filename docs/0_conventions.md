---
layout: default
title: 0. Conventions and useful identities
---


## Form field


$$
\begin{equation}
A = \frac{1}{p!}A_{\mu_1...\mu_p} dx^{\mu_1}...dx^{\mu_p},\quad
A_{\mu_1...\mu_p} = A_{[\mu_1...\mu_p]},\quad
|A|^2 = \frac{1}{p!}A_{\mu_1...\mu_p} A^{\mu_1...\mu_p}=\frac{1}{p!}A^2
\end{equation}
$$


$$
\begin{equation}
(d A)_{\mu_1...\mu_{p+1}} = (p+1)\partial_{[\mu_1} A_{\mu_2...\mu_{p+1}]}
\end{equation}
$$


$$
\begin{equation}
\quad
(A^{(p)} \wedge B^{(q)})_{\mu_1...\mu_p \nu_1 ... \nu_q}
=\frac{(p+q)!}{p!q!}A^{(p)}_{[\mu_1...\mu_p}  B^{(q)}_{\nu_1 ... \nu_q]}
\end{equation}
$$

We will use $\epsilon_{\mu_1 ... \mu_p}$ to denote the Levi-Civita tensor and $\varepsilon_{\mu_1 ... \mu_p}$ to denote the flat-space fully antisymmetrized symbol.


$$
\begin{equation}
\epsilon_{\mu_1...\mu_D} = \sqrt{|g|}\varepsilon_{\mu_1...\mu_D},\quad
\varepsilon_{0,1,...,(D-2),(D-1)} = 1,
\end{equation}
$$


$$
\begin{equation}
\epsilon_{\mu_1...\mu_p \lambda_1 ...\lambda_{D-p}}
\epsilon^{\nu_1...\nu_p \lambda_1 ...\lambda_{D-p}}
=(-1)^{[t]} p!(D-p)!\delta^{\nu_1...\nu_p}_{\mu_1...\mu_p}
\end{equation}
$$


$$
\begin{equation}
(*A)_{\mu_1...\mu_{D-p}} = \frac{1}{p!}
\epsilon_{\mu_1...\mu_{D-p}}{}^{\nu_1...\nu_p} A_{\nu_1...\nu_p},
\end{equation}
$$

## Useful Forms Identities

One can then show that


$$
\begin{equation}
*(*A) = (-1)^{[t]}(-1)^{p(D-p)}A,
\end{equation}
$$


$$
\begin{equation}
(*F)\wedge *(*F)=(-1)^{[t]}F\wedge *F,\quad
|*F|^2 = (-1)^{[t]}|F|^2
\end{equation}
$$

as well as


$$
\begin{equation}
\frac{1}{(D-p-1)!}(*F)_{\mu ...}(*F)_\nu{}^{...}
=(-1)^{[t]}\left[
    \frac{1}{p!}F^2 g_{\mu\nu} - \frac{1}{(p-1)!}F_{\mu...}F_\nu{}^{...}
\right]
\end{equation}
$$

Derivation:


$$
\begin{equation}
\begin{aligned}
\frac{1}{(D-p-1)!}(*F)_{\mu ...}(*F)_\nu{}^{...}
&=\frac{1}{(D-p-1)!}\frac{1}{p!p!}
\epsilon_{\mu \mu_1...\mu_p \lambda_1...\lambda_{D-p-1}}
\epsilon_{\nu} {}^{\nu_1...\nu_p} {}^{\lambda_1...\lambda_{D-p-1}}
F^{\mu_1...\mu_p}F_{\nu_1...\nu_p}\\
&=\frac{p+1}{p!}(-1)^{[t]} g_{\nu\rho} \delta^{[\rho \nu_1...\nu_p]}_{[\mu \mu_1...\mu_p]}F^{\mu_1...\mu_p}F_{\nu_1...\nu_p}
\end{aligned}
\end{equation}
$$


$$
\begin{equation}
\delta^{[\rho \nu_1...\nu_p]}_{[\mu \mu_1...\mu_p]}
\end{equation}
$$

has in total $(p+1)!(p+1)!$ terms and is thus normalized by such number. Within those terms, $(p+1)p!p!$ of which give 


$$
\begin{equation}
\delta^{\rho}_{\mu} \delta^{\nu_1...\nu_p}_{\mu_1...\mu_p}
\end{equation}
$$

equivalent and $(p+1)(p)p!p!$ give 


$$
\begin{equation}
-\delta^{\rho}_{\mu_1}\delta^{\nu_1...\nu_p}_{ \mu ...\mu_p}
\end{equation}
$$

equivalent, so we find


$$
\begin{equation}
\frac{1}{(D-p-1)!}(*F)_{\mu ...}(*F)_\nu{}^{...}
=(-1)^{[t]}\left[\frac{1}{p!}F^2 g_{\mu\nu} - \frac{1}{(p-1)!}F_{\mu...}F_{\nu}{}^{...}\right].
\end{equation}
$$

## Useful gravity identities

$$
\begin{equation}
R(\Lambda^2 g_{mn})
=\Lambda^{-2}[R-2(d-1)\nabla^2\ln\Lambda - (d-2)(d-1)(\nabla \ln \Lambda)^2]
\end{equation}
$$