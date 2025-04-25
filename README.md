# **ParametricDAQP.jl**

**ParametricDAQP.jl** solves multi-parametric quadratic programs of the form

$$
\begin{align}
\min_{z} &  ~\frac{1}{2}z^{T}Hz+(f+F \theta)^{T}z \\
\text{s.t.} & ~A z \leq b + B \theta \\
& ~\theta \in \Theta
\end{align}
$$

where $H \succ 0$ and $\Theta \triangleq \lbrace l \leq \theta \leq u : A_{\theta} \theta \leq b_{\theta}\rbrace$.

The solution $z^*(\theta)$ is a piecewise-affine function over a polyhedral partition.

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
## Citation
If you use the package in your work, consider citing the following paper

```
@inproceedings{arnstrom2024pdaqp,
  author={Arnström, Daniel and Axehill, Daniel},
  booktitle={2024 IEEE 63rd Conference on Decision and Control (CDC)}, 
  title={A High-Performant Multi-Parametric Quadratic Programming Solver}, 
  year={2024},
  volume={},
  number={},
  pages={303-308},
}
```

