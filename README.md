# **ParametricDAQP.jl**

**ParametricDAQP.jl** solves multi-parametric quadratic programs of the form

$$
\begin{align}
\min_{x} &  ~\frac{1}{2}x^{T}Hx+(f+F \theta)^{T}x \\
\text{s.t.} & ~A x \leq b + B \theta \\
& ~\theta \in \Theta
\end{align}
$$

where $H \succ 0$ and $\Theta \triangleq \lbrace l \leq \theta \leq u : A_{\theta} \theta \leq b_{\theta}\rbrace$.

The solution $x^*(\theta)$ is a piecewise-affine function over a polyhedral partition.

## Example
The following code solves the mpQP in Section 7.1 in Bemporad et al. 2002
```julia
using ParametricDAQP

H =  [1.5064 0.4838; 0.4838 1.5258];
f = zeros(2,1);
F = [9.6652 5.2115; 7.0732 -7.0879];
A = [1.0 0; -1 0; 0 1; 0 -1];
b = 2*ones(4);
B = zeros(4,2);
mpQP = (H=H,f=f,F=F,A=A,b=b,B=B);

# Setup parameter region of interest
ub,lb  = 1.5*ones(2), -1.5*ones(2);
Θ = (ub=ub,lb=lb);

# Solve mpQP over desired region
sol,info = ParametricDAQP.mpsolve(mpQP,Θ);
```
