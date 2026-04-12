# Simplifying $t8t8RRDPbDP$
One realizes that
$$
\begin{equation}
\begin{aligned}
t_8 t_8 R^2 |DS^2|
&=\frac{1}{4}{}R^{m n r s} R_{m n r s} S^{t u v w} \bar{S}{}_{t u v w}
+R^{m n r s} R_{m n}{}^{t u} S_{r s}{}^{v w} \bar{S}{}_{t u v w}\\
&-4{}R^{m n r s} R_{m n r}{}^{t} S_{s}{}^{u v w} \bar{S}{}_{t u v w}
-2{}R^{m n r s} R_{m n}{}^{t u} S_{r t}{}^{v w} \bar{S}{}_{s u v w}\\
&+\frac{1}{2}{}R^{m n r s} R^{t u v w} S_{m n r s} \bar{S}{}_{t u v w}
+\frac{1}{2}{}R^{m n r s} R^{t u v w} S_{m n t u} \bar{S}{}_{r s v w}\\
&-4{}R^{m n r s} R_{m}{}^{t u v} S_{n}{}^{w}{}_{u v} \bar{S}{}_{r s t w}
-4{}R^{m n r s} R_{m}{}^{t u v} S_{n}{}^{w}{}_{r s} \bar{S}{}_{t w u v}\\
&-4{}R^{m n r s} R^{t u v w} S_{m n r t} \bar{S}{}_{s u v w}
+8{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} S_{n}{}^{v}{}_{s}{}^{w} \bar{S}{}_{t v u w}\\
&+8{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} S_{n}{}^{v}{}_{u}{}^{w} \bar{S}{}_{s w t v}
+16{}R^{m n r s} R_{m}{}^{t u v} S_{n}{}^{w}{}_{r u} \bar{S}{}_{s v t w}\\
&+4{}R^{m n r s} R^{t u v w} S_{m t r v} \bar{S}{}_{n u s w}
\end{aligned}
\end{equation}
$$

Letting $S_{mnrs} = DP_{mnrs}$, Using Cadabra, we arrive at, here we are using $K_{mn}$ to denote $D_m P_n$, which is symmetric in $m,n$, same for $\bar{K}_{mn}$.

$$
\begin{equation}
\begin{aligned}
t_8 t_8 R^2 |DP^2|
&=\frac{1}{4}{}R^{m n r s} R_{m n r s} K^{t u} \bar{K}{}_{t u}
+\frac{1}{16}{}R^{m n r s} R_{m n r s} K^{t}{}_{t} \bar{K}{}^{u}{}_{u}\\
&-2{}R^{m n r s} R_{m n r}{}^{t} K_{s}{}^{u} \bar{K}{}_{t u}
-R^{m n r s} R_{m n}{}^{t u} K_{r t} \bar{K}{}_{s u}\\
&-\frac{1}{4}{}R^{m n r s} R_{m n r}{}^{t} K^{u}{}_{u} \bar{K}{}_{s t}
-\frac{1}{4}{}R^{m n r s} R_{m n r}{}^{t} K_{s t} \bar{K}{}^{u}{}_{u}\\
&+R^{m n}{}_{m}{}^{r} R^{s t}{}_{s}{}^{u} K_{n r} \bar{K}{}_{t u}
+6{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n t} \bar{K}{}_{s u}\\
&-2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s t u} K_{r t} \bar{K}{}_{s u}
+2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s t u} K_{s t} \bar{K}{}_{r u}\\
&+2{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n s} \bar{K}{}_{t u}
-R^{m n}{}_{m}{}^{r} R^{s t}{}_{s}{}^{u} K_{n t} \bar{K}{}_{r u}\\
&+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r s} K^{t u} \bar{K}{}_{t u}
-2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K_{s}{}^{u} \bar{K}{}_{t u}\\
&+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K^{u}{}_{u} \bar{K}{}_{s t}
+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K_{s t} \bar{K}{}^{u}{}_{u}\\
&+\frac{1}{2}{}R^{m n}{}_{m n} R^{r s t u} K_{r t} \bar{K}{}_{s u}
\end{aligned}
\end{equation}
$$

Sorting

$$
\begin{equation}
\begin{aligned}
t_8 t_8 R^2 |DP^2|
%--- 4 free indices on K \bar{K} ---
&=-{}R^{m n r s} R_{m n}{}^{t u} K_{r t} \bar{K}{}_{s u}
+R^{m n}{}_{m}{}^{r} R^{s t}{}_{s}{}^{u} K_{n r} \bar{K}{}_{t u}\\
&+6{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n t} \bar{K}{}_{s u}
-2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s t u} K_{r t} \bar{K}{}_{s u}\\
&+2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s t u} K_{s t} \bar{K}{}_{r u}
+\cancel{2{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n s} \bar{K}{}_{t u}}\\
&-R^{m n}{}_{m}{}^{r} R^{s t}{}_{s}{}^{u} K_{n t} \bar{K}{}_{r u}
+\frac{1}{2}{}R^{m n}{}_{m n} R^{r s t u} K_{r t} \bar{K}{}_{s u}\\
%--- 2 free indices on K \bar{K} ---
&-2{}R^{m n r s} R_{m n r}{}^{t} K_{s}{}^{u} \bar{K}{}_{t u}
-\frac{1}{4}{}R^{m n r s} R_{m n r}{}^{t} K^{u}{}_{u} \bar{K}{}_{s t}\\
&-\frac{1}{4}{}R^{m n r s} R_{m n r}{}^{t} K_{s t} \bar{K}{}^{u}{}_{u}
-2{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K_{s}{}^{u} \bar{K}{}_{t u}\\
&+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K^{u}{}_{u} \bar{K}{}_{s t}
+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r}{}^{t} K_{s t} \bar{K}{}^{u}{}_{u}\\
%--- 0 free indices on K \bar{K} ---
&+\frac{1}{4}{}R^{m n r s} R_{m n r s} K^{t u} \bar{K}{}_{t u}
+\frac{1}{16}{}R^{m n r s} R_{m n r s} K^{t}{}_{t} \bar{K}{}^{u}{}_{u}\\
&+\frac{1}{2}{}R^{m n}{}_{m}{}^{r} R_{n}{}^{s}{}_{r s} K^{t u} \bar{K}{}_{t u}
\end{aligned}
\end{equation}
$$

One can show that define $K^{mnrs} \equiv D^m\bar{P}^n D^r P^s$.
For any contraction with a real tensor $T$, one has the effective symmetry

$$
\begin{equation}
    K^{mnrs} \cong K^{rsmn}.
\end{equation}
$$

Proof is done by conjugation.
It follows that

$$
\begin{equation}
    K^{mnrs}R_{m}{}^{t}{}_{n}{}^{u}{}R_{rtsu} = 0.
\end{equation}
$$
Apply the first Bianchi identity $R_{r[tsu]}=0$ to the second Riemann factor.

Furthermore,

$$
\begin{equation}
    K^{mnrs}{}R_{mr}{}^{tu}{}R_{nstu}
    = 2{}K^{mnrs}{}R_{m}{}^{t}{}_{r}{}^{u}{}R_{nstu}.
\end{equation}
$$

Proof by First Bianchi identity.

From the properties we have derived, we know that $2{}R^{m n r s} R_{m}{}^{t}{}_{r}{}^{u} K_{n s} \bar{K}{}_{t u} = 0$.
