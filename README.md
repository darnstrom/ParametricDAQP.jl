## **ParametricDAQP.jl**

**ParametricDAQP.jl** solves multi-parametric quadratic programs of the form

$$
\begin{align}
\min_{x} &  ~\frac{1}{2}x^{T}Hx+(f+f_{\theta} \theta)^{T}x \\
\text{s.t.} & ~A x = b + W \theta \\
& ~\theta \in \Theta
\end{align}
$$

where $H \succ 0$ and $\Theta \triangleq \lbrace l \leq \theta \leq u : A_{\theta} \theta \leq b_{\theta}\rbrace$.

The solution $x^*(\theta)$ is a piecewise-affine function over a polyhedral partition.
