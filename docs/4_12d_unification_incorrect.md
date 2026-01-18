---
layout: default
title: 4. 12d Unification with g_s
---


# 12d Unification with g_s

There are several problems that a 12d unification needs to address
- self-duality of $F_5$
- volume of torus
- redundant doFs from higher-dim forms

It would be a miracle if in 12d or in the 12d to 10d reduction ansatz, one is able to address single or multiple of these problems.

## Putting $g_S$ back

In 2512.10746 we had proposed a 12d uplift of the various form-fields when $g_S = 1$.
We note that with $g_S$ restored, as a generic constant, the terms in the type IIB action may be grouped by their couplings to $g_S$ as 
(this form of the action can be obtained by taking e.g. the type IIB action )

$$
\begin{equation}
\begin{aligned}
2\kappa_{10}^2g_S^2S_{IIB}
&=\int d^{10}x\sqrt{-g}\left[
R-\frac{1}{2}\frac{\partial\tau\partial\bar\tau}{\tau_2^2}
\right]\\
&-\frac{g_S}{2}\int d^{10}x\sqrt{-g}\left[
e^{-\Phi}|H_3|^2+e^\Phi|\tilde{F}_3|^2
\right]\\
&-g_S^2\left[\frac{1}{4}\int d^{10}x\sqrt{-g}
|\tilde F_5|^2
+\frac{1}{2}\int C_4\wedge H_3\wedge F_3\right].
\end{aligned}
\end{equation}
$$

The gravity plus axio-dilaton action only depend on the 10d metric and the shape structure of the torus, hence they don't feel the volume of the torus. They also don't couple to $g_S$ (with $\kappa_{10}^2 g_S^2$ factored out).
The SL(2, R) doublet fields should come from wrapping on one-cycles of the torus, hence depend on the volume of the torus, they couple to $g_S$.
The self-dual 5-form may come from wrapping a 7-form completely on a torus, they couple to $g_S^2$.
Suppose the volume of the torus is related to $g_S$, then it is natural to ask whether there is a consistent way to perform the 12d to 10d reduction such that the various $g_S$ couplings arise naturally from wrapping on the torus differently. 
We now explore this possibility.

Let the 12d coordinates be parameterized by $(x^m,u,v)$, with 

$$
\begin{equation}
(u,v)\sim(u+2\pi l,v)\sim (u,v+2\pi l).
\end{equation}
$$ 

Here $l$ has dimension of (length) and is the length scale of the torus.
Let $M_{ij}$ be the $2\times 2$ metric on the torus.
We consider the 12d metric embedding given by

$$
\begin{equation}
\mathcal{G}_{MN}=\begin{pmatrix}
g_{mn} & 0 \\
0 &M_{ij}
\end{pmatrix},\quad
M_{ij}
=\frac{\sqrt{M}}{\tau_2}\begin{pmatrix}
1& \tau_1\\
\tau_1 & \tau_1^2+\tau_2^2
\end{pmatrix}
\end{equation}
$$

where $\det( M_{ij}) = M$ is some dimensionless constant parameterizing the volume:

$$
\begin{equation}
Vol(T_2)=\int_{T_2}\epsilon_2 = \sqrt{M}\int_0^{2\pi l}dudv = 4\pi^2 \sqrt{M}l^2.
\end{equation}
$$

## Axio-dilaton sector
We note that the Ricci scalar evaluated on $\mathcal{G}_{MN}$ is independent of $\sqrt{M}$:

$$
\begin{equation}
\mathcal{R} = R - \frac{\partial\tau\partial\bar\tau}{2\tau_2^2}.
\end{equation}
$$

The axio-dilaton action may be written as

$$
\begin{equation}
\frac{1}{2\kappa_{10}^2g_S^2}\int d^{10}x\sqrt{-g}\left(
R - \frac{\partial\tau \partial\bar\tau}{2\tau_2^2}
\right)
= \underbrace{\left(\frac{1}{2\kappa_{10}^2g_S^2}\frac{M^{-1/2}}{4\pi^2 l^2}\right)}_{2\kappa_{12}^2}
\underbrace{\left(\int_{T_2}^{2\pi l}dudv\int d^{10}x\sqrt{-\mathcal{G}}\mathcal{R}\right)}_{"12d"}.
\end{equation}
$$

Where we have separated 10d parameters and torus parameters such as $\kappa_{10}, g_S, l, M$ and "12d" terms. 
We are hoping that the various 10d and torus parameters reorganize themselves into a 12d parameter $\kappa_{12}$.


## NSNS and RR 3-form sector

The 3-form field strengths form an SL(2, R) doublet. 
To unify them in 12d we define a 12d 4-form field strength with exactly one leg on the torus:

$$
\begin{equation}
\mathcal{F}_4=d\mathcal C_3 = M^{d_4}\cdot \left(H_3\wedge du + F_3\wedge dv\right),
\quad
\mathcal C_3 \equiv M^{d_4}\cdot \left(B_2\wedge du + C_2\wedge dv\right).
\end{equation}
$$

Since we are keeping things general at the moment, we are allowing scaling of the form field by a constant proportional to the volume of torus in the 12 to 10 form field reduction ansatz, given by $M^{d_4}$.


Using the 12d metric ansatz, the 10d 3-form field strengths and their axio-dilaton couplings follow from contracting $\mathcal F_4$:

$$
\begin{equation}
\begin{aligned}
|\mathcal F_4|\bigg|_{\mathcal{G}_{MN}}
&=M^{2d_4}\left(M^{uu}|H_3|^2+2M^{uv}|F_3\cdot H_3|+M^{vv}|F_3|^2\right)\\
&=M^{2d_4-1/2}\left(
e^{-\Phi}|H_3|^2 + e^\Phi|F_3-CH_3|^2
\right)\bigg|_{g_{mn}}.
\end{aligned}
\end{equation}
$$

Hence we find a 12d interpretation of the NSNS and RR 3-form sector:

$$
\begin{equation}
\frac{1}{2\kappa_{10}^2 g_S^2}\int d^{10}x \sqrt{-g}\left(-\frac{g_S}{2}\right)
\left[e^{-\Phi}|H_3|^2+e^\Phi|\tilde{F}_3|^2\right]
=\frac{1}{2\kappa_{10}^2 g_S^2}\frac{g_SM^{-2d_4}}{4\pi^2l^2}\int_{T_2} dudv \int d^{10}x\sqrt{-\mathcal{G}}\left[-\frac{1}{2}|\mathcal{F}_4|^2\right].
\end{equation}
$$

## Self-dual 5-form

Finding a 12d interpretation for the 5-form sector is trickier. 
One observes that the composite, SL(2, R) singlet 5-form 

$$
\begin{equation}
\tilde F_5 = F_5 - \frac12 C_2 \wedge H_3 + \frac{1}{2}B_2\wedge F_3    
\end{equation}
$$

can not be sensibly constructed in 12d, due to the lack of a pair of form field potential and strength with ranks that sum to 5 in 12d.
However, the 10d self-dual 5-form field strength admits two possible uplifts to 12d, either as a 5-form or 7-form:

$$
\begin{equation}
\mathcal{C}_4 = M^{d_5} C_4,\quad \mathcal{F}_5 = M^{d_5} F_5,
\end{equation}
$$

$$
\begin{equation}
\mathcal{F}_7=M^{d_7}\cdot F_5\wedge \epsilon_2 = M^{d_7+1/2}\cdot F_5\wedge du\wedge dv.
\end{equation}
$$

Then note that

$$
\begin{equation}
M^{-2d_7}|\mathcal{F}_7|^2\bigg|_{\mathcal{G}_{MN}}
= M^{-2d_5}|\mathcal{F}_5|^2\bigg|_{\mathcal{G}_{MN}} 
= |F_5|^2\bigg|_{g_{mn}},
\end{equation}
$$

which allows us to define a 12d composite 7-form (consistency here requires $2d_4 = d_7+1/2$ )

$$
\begin{equation}
\tilde{\mathcal F}_7 \equiv \mathcal F_7 + \frac{1}{2}\mathcal{C}_3 \wedge \mathcal F_4
\end{equation}
$$

that supplies the type IIB composite 5-form contribution upon contraction in 12d:

$$
\begin{equation}
M^{-2d_7}|\tilde{\mathcal{F}}_7|^2\bigg|_{\mathcal{G}_{MN}} = |\tilde F_5|^2\bigg|_{g_{mn}}.
\end{equation}
$$

Hence

$$
\begin{equation}
\frac{1}{2\kappa_{10}^2g_S^2}\left(-\frac{g_S^2}{4}\right)
\int d^{10}x\sqrt{-g}|\tilde{F}_5|^2
=\frac{1}{2\kappa_{10}^2g_S^2}\frac{g_S^2M^{-2d_7-1/2}}{4\pi^2l^2}\int_{T_2} dudv\int d^{10}x\sqrt{-\mathcal{G}}\left[-\frac{1}{4}|\tilde{\mathcal{F}}_7|^2\right]_{\mathcal{G}_{MN}}
\end{equation}
$$


The 10d self-duality condition on $F_5$ may then be written as a 12d Hodge duality

$$
\begin{equation}
M^{d_7-d_5}\mathcal F_5 = *_{12}\mathcal F_7\quad
\Leftrightarrow\quad
F_5 = *_{10}F_5.
\end{equation}
$$

It is very tempting to impose $d_5 = d_7$, since that would relate the 5- and 7-form algebraically by the 12d Hodge duality.

## Chern-Simons term

The 10d Chern-Simons term can also be obtained using the 12d potentials and their corresponding field strengths:

$$
\begin{equation}
C_4\wedge H_3\wedge F_3\wedge du\wedge dv = \frac{1}{2}M^{-d_5-2d_4}\mathcal{C}_4\wedge\mathcal{F}_4\wedge \mathcal{F}_4
\end{equation}
$$

which implies

$$
\begin{equation}
\frac{1}{2\kappa_{10}^2 g_S^2}(-g_S^2)\frac{1}{2}\int C_4\wedge H_3\wedge F_3
=\frac{1}{2\kappa_{10}^2g_S^2}\frac{g_S^2M^{-d_5-2d_4}}{4\pi^2l^2}\left(-\frac{1}{4}\right)\int_{T_2\times \mathbb{R}^{1,9}}\mathcal{C}_4\wedge \mathcal{F}_4\wedge \mathcal{F}_4
\end{equation}
$$

## A 12d Unificaton

We expect the 12d unification to be given in the form of an integration over $T_2\times \mathbb{R}^{1,9}$ multiplying 10d and torus quantities $g_S, l_S, \sqrt{M}$, which reorganizes themselves into some 12d quantity $\kappa_{12}^2$ given solely by a 12d length scale $l_F$.
That is, $\kappa_{12}^2 = l_F^{10}$ ignoring factors of $2\pi$.

For the lack of a third scalar in type IIB the volume  (actually area) of the torus shall be fixed.
Drawing from our experience on the type IIA side, the constant volume of the torus may be related to $g_S$.
Then the various sectors discussed above, which arise from different ways of wrapping 12d fields on the torus, will pick up nontrivial factors of $g_S$.
From our 12d to 10d dimensional reduction ansatz we also allowed scalings of the form fields by constants proportional to the volume of the torus $M^{d_4}, M^{d_5}, M^{d_7}$, where consistency requries $2d_4=d_7+1/2$ and $d_5 = d_7$. 
They will also introduce factors of $g_S$ as $\sqrt{M} \sim Vol(T_2)$ is related to $g_S$.
It is thus worthwhile to investigate, with when all these sources of $g_S$ are taken together, is it possible to obtain the specific $g_S$ couplings of the various form fields in the type IIB action?

This is indeed is the case, with

$$
\begin{equation}
\begin{aligned}
S_{IIB}=
S_{\text{``12"}}
= \frac{1}{2\kappa_{12}^2}\int_{T_2}dudv\int d^{10}x
\sqrt{-\mathcal G}
\left(
\mathcal{R} - \frac{1}{2}|\mathcal F_4|^2-\frac{1}{4}|\tilde{\mathcal{F}}_7|^2
\right)-\frac{1}{4\kappa_{12}^2}\int_{T_2\times \mathbb{R}^{1,9}}\mathcal C_4\wedge \mathcal F_4\wedge \mathcal F_4,
\end{aligned}
\end{equation}
$$


$$
\begin{equation}
\frac{1}{2\kappa_{12}^2}=\frac{1}{2\kappa_{10}^2g_S^2}\frac{1}{4\pi^2\sqrt{M} l^2}=\frac{1}{2\kappa_{10}^2 g_S^2 Vol(T_2)} ,\quad
M^{d_7}=g_S, \quad
d_5 = d_7, \quad
2d_4 = d_7 + \frac{1}{2}.
\end{equation}
$$

As a quick reminder the integration on $T_2$ is between 0 and $2\pi l$ whereas $l$ is understood as the length unit on the torus.

Now what is $l$? 
It does not matter, because we are free to measure length on the torus with whatever length unit we'd like, the numeric, dimensionless measurements just have to change accordingly.
For example, the dimensionless volume $\sqrt{M}$ depends on what we are using for $l$:

$$
\begin{equation}
\sqrt{M}l^2 = \sqrt{M_S}l_S^2 = \sqrt{M_F}l_F^2.
\end{equation}
$$

Suppose we choose $l = l_F$, then we obtain the relation

$$
\begin{equation}
l_F^{8} = \sqrt{M_F}g_S^2 l_S^8.
\end{equation}
$$

If one believes that $l = l_F$, i.e. that $\sqrt{M_F} = 1$, they would find

$$
\begin{equation}
l_F^4 = g_Sl_S^4,\quad
\sqrt{M_S} = g_S^{1/2},\quad
Vol(T_2) = 4\pi^2 g_S^{1/2}l_S^2.
\end{equation}
$$