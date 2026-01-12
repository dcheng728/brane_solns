## Setup, Conventions, and EQMs
Much of this follows [^Stelle98].

In $D$ dimensions, for the Einstein-Hilbert action coupled to various $n$-forms with no dilaton,

$$
S = \int d^D x \sqrt{-g}\left[R - \sum_n\frac{1}{2}\frac{1}{n!}F_n^2\right]
\qquad(1)
$$

The Einstein equations are

$$
\begin{aligned}
R_{MN} 
= S_{MN} 
&= \sum_n \frac{1}{2(n-1)!}
\left[
(F_n)_{M...}(F_n)_N{}^{...}-\frac{n-1}{n(D-2)}(F_n)^2g_{MN}  
\right]\\
&=\frac{1}{2}\sum_n \left[
    \frac{1}{(n-1)!}F_{M...}F_N^{...}
    -\frac{n-1}{D-2}\left(\frac{1}{n!}F^2\right)g_{MN}
\right],
\end{aligned}
$$

$$
dF = 0.
$$

We look for solutions that are $d$-dimensional Minkowski space transverse to $(D-d)$-dimensional Euclidean space.
So our coordinates are split into $x^M=(x^\mu,y^m)$, where $x^\mu$ for $\mu = 0,1,...,d-1$ are coordinates on the Minkowski space, and $y^m$ for $m = d,d+1,...,D-1$ are embedding coordinates on the $D-d$ dimensional Euclidean space.



### Form field conventions

$$
A^{(p)} = \frac{1}{p!}A^{(p)}_{\mu_1...\mu_p} dx^{\mu_1}...dx^{\mu_p},\quad
A^{(p)}_{\mu_1...\mu_p} = A^{(p)}_{[\mu_1...\mu_p]},
$$

$$
(d A^{(p)})_{\mu_1...\mu_{p+1}} = (p+1)\partial_{[\mu_1} A^{(p)}_{\mu_2...\mu_{p+1}]}
$$

$$
\quad
(A^{(p)} \wedge B^{(q)})_{\mu_1...\mu_p \nu_1 ... \nu_q}
=\frac{(p+q)!}{p!q!}A^{(p)}_{[\mu_1...\mu_p}  B^{(q)}_{\nu_1 ... \nu_q]}
$$

We will use $\epsilon_{\mu_1 ... \mu_p}$ to denote the Levi-Civita tensor and $\varepsilon_{\mu_1 ... \mu_p}$ to denote the flat-space fully antisymmetrized symbol.

$$
\epsilon_{\mu_1...\mu_D} = \sqrt{|g|}\varepsilon_{\mu_1...\mu_D},\quad
\varepsilon_{01...(D-2)(D-1)} = 1,
$$

$$
(*A)_{\mu_1...\mu_{D-p}} = \frac{1}{p!}
\epsilon_{\mu_1...\mu_{D-p}}{}^{\nu_1...\nu_p} A^{(p)}_{\nu_1...\nu_p},
$$

$$
*(*A) = (-1)^{[t]}(-1)^{p(D-p)}A.
$$

## Brane-solution ansatz and its implications
We usually solve with the ansatz for a $(p=d-1)$-brane:

$$
ds_D^2 = H^{a}(r)\eta_{\mu\nu}dx^\mu dx^\nu + H^b(r) \delta_{mn}dy^mdy^n,\quad
r\equiv \sqrt{y^my^m}.
$$

where $H$ is a harmonic function in $(D-d)$-dimensions, usually taken to be

$$
H = 1 + \frac{Q}{r^{\tilde{d}}},\quad \tilde{d} = D-d-2
$$

Then the Ricci tensor is found to be

$$
R_{\mu\nu}
=-\eta_{\mu\nu}(H')^2H^{a-b-2}\left[\frac{a(ad+b\tilde d-2)}{4}\right]
$$

<!-- $$
\mathcal{B} = \frac{(H')^2}{4H^2}b(ad+b\tilde{d}-2)+\frac{H'}{rH}\frac{(ad+b\tilde{d})}{2}
$$

$$
\mathcal{C} = \frac{(H')^2}{4H^2}[ad(a-b)-(ad+b\tilde{d})(b+2)]
+\frac{H'}{rH}\frac{(ad+b\tilde{d})}{2}(-2-\tilde{d})
$$ -->

$$
\begin{aligned}
R_{mn} 
&= \frac{(H')^2}{4H^2}\left[-b(ad+b\tilde{d}-2)\delta_{mn}-\frac{y^m y^n}{r^2}\left[ad(a-b)-(ad+b\tilde{d})(b+2)\right]\right]\\
&+\frac{H'}{rH}\frac{(ad+b\tilde{d})}{2}\left[-\delta_{mn} + \frac{y^m y^n}{r^2}(2+\tilde{d})\right]
\end{aligned}
$$

Usually one adopts the ansatz that $A_{n-1} \propto H$, hence $F_n \propto H'$.
Then $S_{MN}$ will be exactly quadratic in $H'$ up to some powers of $H$: 

$$S_{MN}\sim H^{(...)} (H')^2.$$

which means in order to satisfy the Einstein equation, $R_{MN}$ should be exactly quadratic in $H'$.
This demands

$$
b=-ad/\tilde{d},
$$

and we can write down a simplified Ricci tensor

$$
R_{\mu\nu} = \frac{a}{2}(H')^2H^{a-b-2} \eta_{\mu\nu}
=\frac{a}{2}(H')^2H^{\frac{ad+(a-2)(D-d-2)}{(D-d-2)}} \eta_{\mu\nu},
$$

$$
R_{mn} =\frac{(H')^2}{4H^2}\left[2b\delta_{mn}-ad(a-b)\frac{y^my^n}{r^2}\right]
=-\frac{ad}{2(D-d-2)}\frac{(H')^2}{H^2}\left[
\delta_{mn}+\frac{a(D-2)}{2}\frac{y^m y^n}{r^2}
\right]
$$

## Electric ansatz
For a $F_{d+1}$ form field, the electric ansatz can be given by

$$
(F_{d+1})_{m \mu_0 \mu_1...\mu_{d-1}} = c\varepsilon_{\mu_0 \mu_1...\mu_{d-1}}H^{-1 + \frac{ad}{2}} \partial_m H
$$

which satisfies $d F_{d+1} = 0$. The exponent $-1 + \frac{ad}{2}$ is required by matching the scaling in $H$ on both sides of $R_{mn} = S_{mn}$.
One can then solve to find that

$$
a = -\frac{2}{d},\quad
b = \frac{2}{\tilde{d}},\quad
c^2 = \frac{2(D-2)}{d\tilde{d}}.
$$

There are then no free parameters, and we end up with the usual electric brane solutions as described in [^Stelle98].

## Magnetic ansatz
For a $F_{d+1}$ form field, the electric ansatz can be given by




## 12d to 10d reduction ansatz
What should be the reduction ansatz from 12 to 10?
In 12d one should not know about the torus, hence know nothing about 10d.
Only once the torus is introduced the 10d theory is obtained via 12 = 10 + 2.
The reduction ansatz should only break symmetry as 12 = 10 + 2.

$$
ds_{12}^2 = H^f ds_2^2 + H^g ds_{10}^2
$$

$$
ds_{12}^2 = H^f ds_2^2 + H^{g-\frac{1}{2}}ds_{1,3}^2 + H^{g+\frac{1}{2}}ds_6^2
$$

Then either (1) the reduction came from 12d 5-brane, $f=g-\frac{1}{2}$ or (2) the reduction came from 12d 3-brane $f=g+\frac{1}{2}$.

### 12d 5-brane

$$
ds_{12}^2 = H^f ds_{1,5}^2 + H^{f+1}ds_6^2
$$

Solving using above formula yields

$$
ds_{12}^2 = H^{-\frac{2}{5}}ds_{1,5}^2 +H^{\frac{3}{5}}ds_6^2
$$

### 12d 3-brane

$$
ds_{12}^2 = H^{f-1} ds_{1,3}^2 + H^{f}ds_8^2
$$

Solving using above formula yields

$$
ds_{12}^2 = H^{-\frac{3}{5}}ds_{1,3}^2 +H^{\frac{2}{5}}ds_8^2
$$

This makes sense; try to balance with this.

### Alternative

Following Kelly's paper, you shall find 12d metrics balanced by 5 and 7-forms

$$
ds_{12}^2 = H^{-1/2}ds_{1,3}^2 + H^{2/3}ds_8^2
$$

$$
ds_{12}^2 = H^{-1/3}ds_{1,5}^2 + H^{1/2}ds_6^2
$$

But if you wanted to reduce these to 10d and obtain the D3 solution, the 12 to 10 ansatz is not in 12=2+10 form. So perhaps you should not try to follow this.



[^Stelle98]: K. S. Stelle, 
*BPS branes in supergravity*,
*ICTP Summer School in High-Energy Physics and Cosmology* (1997),
[arXiv:hep-th/9803116](https://inspirehep.net/literature/468167).