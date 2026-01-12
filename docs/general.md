## Setup, Conventions, and EQMs
Much of this is follows [Stelle, 9803116]

In $D$ dimensions, for Einstein-Hilbert action coupled to various n-forms with no dilaton

$$
S = \int d^D x \sqrt{-g}\left[R - \sum_n\frac{1}{2}\frac{1}{n!}F_n^2\right]
$$

The Einstein equation is

$$
R_{MN} 
= S_{MN} 
= \sum_n \frac{1}{2(n-1)!}
\left[
(F_n)_{M...}(F_n)_N{}^{...}-\frac{n-1}{n(d+\tilde{d})}(F_n)^2g_{MN}  
\right]
$$

$$
dF = 0
$$

We look for solutions that are $d$ dimensional Minkowski space transverse to $D-d$ dimensional Euclidean space.
So our coordinates are split into $x^M=(x^\mu,y^m)$, where $x^\mu$ for $\mu = 0,1,...,d-1$ are coordinates on the Minkowski space, and $y^m$ for $m = d,d+1,...,D-1$ are embedding coordinates on the $D-d$ dimensional Euclidean space.

## Brane-solution implications
We usually solve with the ansatz for $(d-1)$-brane

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

$$
\mathcal{B} = \frac{(H')^2}{4H^2}b(ad+b\tilde{d}-2)+\frac{H'}{rH}\frac{(ad+b\tilde{d})}{2}
$$

$$
\mathcal{C} = \frac{(H')^2}{4H^2}[ad(a-b)-(ad+b\tilde{d})(b+2)]
+\frac{H'}{rH}\frac{(ad+b\tilde{d})}{2}(-2-\tilde{d})
$$

$$
\begin{aligned}
R_{mn} 
&= \frac{(H')^2}{4H^2}\left[-b(ad+b\tilde{d}-2)\delta_{mn}-\frac{y^m y^n}{r^2}\left[ad(a-b)-(ad+b\tilde{d})(b+2)\right]\right]\\
&+\frac{H'}{rH}\frac{(ad+b\tilde{d})}{2}\left[-\delta_{mn} + \frac{y^m y^n}{r^2}(2+\tilde{d})\right]
\end{aligned}
$$

We will also use ansatzs that $A_{n-1} \propto H$ hence $F_n \propto H'$, for some harmonic function $H$.
Then it's safe to state that $S_{MN}$ will be exactly quadratic in $H'$: 

$$S_{MN}\sim H^{(...)} (H')^2.$$

which means after substituting the harmonic condition of $H$, in order to satify the Einstein equation $R_{MN}$ should be exactly quadratic in $H'$.
This allows us to obtain

$$
b=-a\left(\frac{d}{D-d-2}\right).
$$

## 12d to 10d
What should be the reduction ansatz from 12 to 10?

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

This makes sense, try to balance with this!

### Alternative

Following Kelly's paper, you shall find 12d metrics balanced by 5 and 7-forms

$$
ds_{12}^2 = H^{-1/2}ds_{1,3}^2 + H^{2/3}ds_8^2
$$

$$
ds_{12}^2 + H^{-1/3}ds_{1,5}^2 + H^{1/2}ds_6^2
$$

But if you wanted to reduce these to 10d and obtain the D3 solution, the 12 to 10 ansatz is not in 12=2+10 form. So perhaps you should not try to follow this.
