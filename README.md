# SVE-simulation
This repository contains simulation code for stochastic Volterra equations (SVE).
We consider the stochastic Volterra equation

$$
 X_t = X_0
    + \int_0^t K(t - s)\,b\bigl(X_s\bigr)\,\mathrm{d}s
    + \int_0^t K(t - s)\,\sigma\bigl(X_s\bigr)\,\mathrm{d}W_s,
$$

where $K(u) = \frac{u^{H-\frac12}}{\Gamma\bigl(H+\tfrac12\bigr)}, \ H\in\bigl(0,\tfrac12\bigr].$ We approximate the solution of the SVE by a numerical scheme based on a discrete sampling grid. The (possibly non-equidistant) time grid is defined recursively as

$$
 \tau_0^n=0,
  \quad
  \tau_{k+1}^n
  =\tau_k^n+\frac{1}{n\theta(\tau_k^n)},
  $$

  where $\theta_t$ is a strictly positive $c\grave{a}dl\grave{a}g$ adapted process and $|1/\theta_t|$ is bounded by a constant $K_0$. 

  Both equidistant and non-equidistant sampling schemes are considered in the numerical implementation. Under suitable regularity assumptions, the numerical approximation is expected to achieve the convergence rate $n^{-H}$. To verify the theoretical convergence rate, we perform numerical simulations using MATLAB. The experiments are conducted for different values of the Hurst parameter H=0.1,0.3,0.4. For simplicity, the drift and diffusion coefficients are chosen as $b=\sigma=1$ and the sampling term $\theta$ as

$$
\theta_t = 1 + 0.5|W_t|.
$$
