# SVE-simulation

This repository contains simulation code for stochastic Volterra equations (SVEs).

We consider the stochastic Volterra equation


$$X_t=X_0+\int_0^t K(t-s)b(X_s)\,ds+\int_0^t K(t-s)\sigma(X_s)\,dW_s,$$

where
$$K(u)=\frac{u^{H-\frac12}}{\Gamma\left(H+\frac12\right)},\qquad H\in\left(0,\frac12\right].$$

We approximate the solution of the SVE by a numerical scheme based on a discrete sampling grid. The (possibly non-equidistant) time grid is defined recursively by

$$\tau_0^n=0,\qquad\tau_{k+1}^n=\tau_k^n+\frac{1}{n\theta(\tau_k^n)},$$

where $\theta_t$ is a strictly positive adapted càdlàg process and
$|1/\theta_t|$ is bounded by a constant $K_0$.

Under the assumptions of the theoretical result, the numerical approximation is predicted to achieve the convergence rate

$$
E_n = O(n^{-H}),
$$

where $E_n$ denotes the approximation error.

To verify this theoretical convergence rate, we perform numerical simulations using MATLAB. The experiments are conducted for different values of the Hurst parameter,

$$
H=0.1,\quad 0.3,\quad 0.4.
$$

For simplicity, the drift and diffusion coefficients are chosen as

$$
b=\sigma=1,
$$

and the sampling process is specified by

$$
\theta_t=1+0.5|W_t|.
$$

This choice gives a strictly positive continuous (and hence càdlàg) adapted process and satisfies

$$
\frac{1}{\theta_t}\leq 1.
$$

The convergence rate is assessed on a log-log scale by plotting the approximation error $E_n$ against $n$ and estimating the slope through linear regression. According to the theoretical result,

$$
\log E_n
\approx
C-H\log n,
$$

so that the theoretical slope is $-H$.
