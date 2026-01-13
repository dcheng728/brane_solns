# A 12d Unification

There are several problems that a 12d unification needs to address
- self-duality of $F_5$
- volume of torus
- redundant doFs from higher-dim forms

It would be a miracle if in 12d or in the 12d to 10d reduction ansatz, one is able to address single or multiple of these problems.

---
We first note that with $g_s$ restored, the terms in the type IIB action may be grouped by their couplings to $g_s$ as

$$
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
$$


Let the 12d coordinates be parameterized by $(x^m,u,v)$, and let $M_{ij}$ be the $2\times 2$ metric on the torus.
One may consider the 12d metric embedding given by

$$
\mathcal{G}_{MN}=\begin{pmatrix}
g_{mn} & 0 \\
0 &M_{ij}
\end{pmatrix},\quad
M_{ij}
=\frac{\sqrt{M}}{\tau_2}\begin{pmatrix}
1& \tau_1\\
\tau_1 & \tau_1^2+\tau_2^2
\end{pmatrix}
$$

where $\det( M_{ij}) = M$ is some constant.

Then the volume of the torus can be obtained by integrating the volume form $\epsilon_2 = \sqrt{M}du\wedge dv$:

$$
Vol(T_2)=\int_{T_2}\epsilon_2 = \sqrt{M}\int_0^{2\pi l}dudv = 4\pi^2 \sqrt{M}l^2,
$$

where $l$ is the length scale of the torus. 
The axio-dilaton action may be written 
In particular,

$$
\frac{1}{2\kappa_{10}^2g_s^2}\int d^{10}x\sqrt{-g}\left(
R - \frac{\partial\tau \partial\bar\tau}{2\tau_2^2}
\right)
= \underbrace{\left(\frac{1}{2\kappa_{10}^2g_s^2}\frac{M^{-1/2}}{4\pi^2 l^2}\right)}_{2\kappa_{12}^2}
\underbrace{\left(\int_{T_2}dudv\int d^{10}x\sqrt{-\mathcal{G}}\mathcal{R}\right)}_{"12d"}.
$$

Where we have separated 10d parameters and torus parameters such as $\kappa_{10}, g_s, l, M$ on one side and "12d" terms on the other. 
In 12d one should not know about $\kappa_{10}$, or the torus at all.
We are hoping that the various 10d and torus parameters reorganize themselves into a 12d parameter $\kappa_{12}$.


The 3-form field strengths form an SL(2, R) doublet. 
To unify them in 12d we define a 12d 4-form field strength with exactly one leg on the torus:

$$
\mathcal{F}_4=d\mathcal C_3 = M^{d_4}\cdot \left(H_3\wedge du + F_3\wedge dv\right),
\quad
\mathcal C_3 \equiv M^{d_4}\cdot \left(B_2\wedge du + C_2\wedge dv\right).
$$

Since we are keeping things general at the moment, we are allowing scaling of the form field by a constant proportional to the volume of torus in the 12 to 10 form field reduction ansatz.


Using the 12d metric ansatz \eqref{12d metric ansatz}, the 10d 3-form field strengths and their axio-dilaton couplings follow from contracting $\mathcal F_4$:

$$
\begin{aligned}
|\mathcal F_4|\bigg|_{\mathcal{G}_{MN}}
&=M^{2d_4}\left(M^{uu}|H_3|^2+2M^{uv}|F_3\cdot H_3|+M^{vv}|F_3|^2\right)\\
&=M^{2d_4-1/2}\left(
e^{-\Phi}|H_3|^2 + e^\Phi|F_3-CH_3|^2
\right)\bigg|_{g_{mn}}.
\end{aligned}
$$

Hence

$$
\frac{1}{2\kappa_{10}^2 g_s^2}\int d^{10}x \sqrt{-g}\left(-\frac{g_s}{2}\right)
\left[e^{-\Phi}|H_3|^2+e^\Phi|\tilde{F}_3|^2\right]
=\frac{1}{2\kappa_{10}^2 g_s^2}\frac{g_sM^{-2d_4}}{4\pi^2l^2}\int_{T_2} dudv \int d^{10}x\sqrt{-\mathcal{G}}\left[-\frac{1}{2}|\mathcal{F}_4|^2\right]
$$

Finding a 12d interpretation for the 5-form sector is trickier. 
One observes that the composite, SL(2, R) singlet 5-form 

$$
\tilde F_5 = F_5 - \frac12 C_2 \wedge H_3 + \frac{1}{2}B_2\wedge F_3    
$$

can not be sensibly constructed in 12d, due to the lack of a pair of form field potential and strength with ranks that sum to 5 in 12d.
However, the 10d self-dual 5-form field strength admits two possible uplifts to 12d:

$$
\mathcal{C}_4 = M^{d_5} C_4,\quad \mathcal{F}_5 = M^{d_5} F_5,
$$

$$
\mathcal{F}_7=M^{d_7}\cdot F_5\wedge \epsilon_2 = M^{d_7+1/2}\cdot F_5\wedge du\wedge dv.
$$

Then note that

$$
M^{-2d_7}|\mathcal{F}_7|^2\bigg|_{\mathcal{G}_{MN}}
= M^{-2d_5}|\mathcal{F}_5|^2\bigg|_{\mathcal{G}_{MN}} 
= |F_5|^2\bigg|_{g_{mn}},
$$

which allows us to define a 12d composite 7-form (suggesting $2d_4 = d_7+1/2$ )

$$
\tilde{\mathcal F}_7 \equiv \mathcal F_7 + \frac{1}{2}\mathcal{C}_3 \wedge \mathcal F_4
$$

that supplies the type IIB composite 5-form contribution upon contraction in 12d:

$$
M^{-2d_7}|\tilde{\mathcal{F}}_7|^2\bigg|_{\mathcal{G}_{MN}} = |\tilde F_5|^2\bigg|_{g_{mn}}.
$$

Hence

$$
\frac{1}{2\kappa_{10}^2g_s^2}\left(-\frac{g_s^2}{4}\right)
\int d^{10}x\sqrt{-g}|\tilde{F}_5|^2
=\frac{1}{2\kappa_{10}^2g_s^2}\frac{g_s^2M^{-2d_7-1/2}}{4\pi^2l^2}\int_{T_2} dudv\int d^{10}x\sqrt{-\mathcal{G}}\left[-\frac{1}{4}|\tilde{\mathcal{F}}_7|^2\right]_{\mathcal{G}_{MN}}
$$


The 10d self-duality condition on $F_5$ may then be written as a 12d Hodge duality

$$
M^{d_7-d_5}\mathcal F_5 = *_{12}\mathcal F_7\quad
\Leftrightarrow\quad
F_5 = *_{10}F_5.
$$

If one wishes to define $\mathcal{F}_7$ this way, then $\mathcal{F}_5$ and $\mathcal{F}_7$ are not independent doFs, they are just two equivalent ways to write down the same underlying doF. 
If one indeed chooses to do that, since $d_5, d_7$ appear in the 12 to 10 reduction ansatz, it is very suggestive that $d_5 = d_7$.

The 10d Chern-Simons term can also be obtained using the 12d potentials and their corresponding field strengths:

$$
C_4\wedge H_3\wedge F_3\wedge du\wedge dv = \frac{1}{2}M^{-d_5-2d_4}\mathcal{C}_4\wedge\mathcal{F}_4\wedge \mathcal{F}_4
$$

which implies

$$
\frac{1}{2\kappa_{10}^2 g_s^2}(-g_s^2)\frac{1}{2}\int C_4\wedge H_3\wedge F_3
=\frac{1}{2\kappa_{10}^2g_s^2}\frac{g_s^2M^{-d_5-2d_4}}{4\pi^2l^2}\left(-\frac{1}{4}\right)\int_{T_2\times \mathbb{R}^{1,9}}\mathcal{C}_4\wedge \mathcal{F}_4\wedge \mathcal{F}_4
$$

We expect the 12d action to be written in the form of integration over $T_2\times \mathbb{R}^{1,9}$ multiplying some parameters that live in 10d which reorganizes into some 12d parameter.
This is indeed possible, one can write down a "12d" action with the correct $g_s$ couplings in 10d

$$
\begin{aligned}
S_{IIB}=
S_{\text{``12"}}
= \frac{1}{2\kappa_{12}^2}\int_{T_2}dudv\int d^{10}x
\sqrt{-\mathcal G}
\left(
\mathcal{R} - \frac{1}{2}|\mathcal F_4|^2-\frac{1}{4}|\tilde{\mathcal{F}}_7|^2
\right)-\frac{1}{4\kappa_{12}^2}\int_{T_2\times \mathbb{R}^{1,9}}\mathcal C_4\wedge \mathcal F_4\wedge \mathcal F_4.
\end{aligned}
$$


with 

$$
\frac{1}{2\kappa_{12}^2}=\frac{1}{2\kappa_{10}^2}\frac{1}{4\pi^2 l^2} \frac{1}{g_s^2 \sqrt{M}}
$$

and $M^{d_7}=g_s$, $d_5 = d_7$, $2d_4 = d_7 + \frac{1}{2}$, and $d_7$ not yet determined (it would be really nice to find a way to constrain it).

As a quick reminder $l$ was defined to be the unit of length of the torus.
A key question now is how should $l$ be?
It should not matter, because after all it's just a choice of in what unit did we measure it in.
For us, there is no other choice but to use $l = l_s$.
To avoid any ambiguity in the consideration of length scales, we will let $\sqrt{M_s}$ be $\sqrt{M}$ measured  with $l_s$, i.e. $vol(T_2) = 4\pi^2\sqrt{M_s}l_s^2$.

In 12d there should be a length scale $l_F$, since there shouldn't be any parameters besides a length scale, we expect (we now ignore constant factors for simplicity)

$$
\frac{1}{\kappa_{12}^2} = \frac{1}{l_F^{10}} = \frac{1}{\kappa_{10}^2 g_s^2 \sqrt{M_s}}
$$

which suggests 

$$
\sqrt{M_s} = \frac{l_F^{10}}{l_s^{10}g_s^2}.
$$


<!-- 


The $\mathcal F_7$ does not arise from an independent degree of freedom, it is related to $\mathcal C_4$.
The action \eqref{12d action} is exactly the Einstein frame type IIB action \eqref{Einstein frame IIB action}.
The 2d integrand is just a repackaging of $1/\kappa_{10}^2$, and it is likely that $\mathcal C_3, \mathcal C_4$ do not furnish independent degrees of freedom, as has been discussed in \cite{Tseytlin:1996ne}. -->
