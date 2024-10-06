# **ParametricDAQP.jl**

**ParametricDAQP.jl** solves multi-parametric quadratic programs of the form

$$
\begin{align}
\min_{x} &  ~\frac{1}{2}x^{T}Hx+(f+f_{\theta} \theta)^{T}x \\
\text{s.t.} & ~A x \leq b + W \theta \\
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
f_theta = [9.6652 5.2115; 7.0732 -7.0879];
A = [1.0 0; -1 0; 0 1; 0 -1];
W = zeros(4,2);
b = 2*ones(4);
senses = zeros(Cint,5); # inequality constraints => sense = 0
mpQP = (H=H,f=f,f_theta=f_theta,A=A,b=b,W=W,senses=senses);

# Setup parameter region of interest
ub,lb  = 1.5*ones(2), -1.5*ones(2);
Θ = (A = zeros(2,0), b=zeros(0), ub=ub,lb=lb);

# Solve mpQP over desired region
opts = ParametricDAQP.EMPCSettings();
F,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
```
