# Fractional-nonlinear-Schrödinger-equation

This notebook shows a simple, scalar implementation of the [finite difference](https://en.wikipedia.org/wiki/Finite_difference_method) for solving the [Nonlinear Schrödinger Equation](https://en.wikipedia.org/wiki/Nonlinear_Schrödinger_equation).

$\frac{\partial A}{\partial z}=-\frac{\alpha}{2}A+i \frac{\beta_2}{2} \frac{\partial^2 A}{\partial t^2}-i \gamma(|A|^2A)$

This nonlinear partial differential equation models how the envelope and phase of light pulse changes when propagating through a single mode optical fiber, when taking power attenuation ($\alpha$), group velocity dispersion ($\beta_2$)and waveguide nonlinearity ($\gamma$) causing self-phase modulation (SPM) into account. A is the slowly varying amplitude of the pulse envelope and t is measured in a frame of reference moving with the pulse at the group velocity $v_g$. The nonlinear Schrödinger equation (for optics) can be derived from the wave equation. However we can choose between two (Fourier transform) notations. The derivation with the negative one can be found in Ursula Keller's [book](https://link.springer.com/book/10.1007/978-3-030-82532-4). I used this, because the scipy library uses the negative convention for the DFT. Depending on the initial width $T_0$ and the peak power $P_0$ of the incident pulse, either dispersive or non linear
effects may dominate along the fiber. It is useful to introduce two length scales, known as the dispersion length $L_D$ and the nonlinear length $L_{NL}$. Let us consider a time scale normalized to the input width $T_0$ as:

$\tau=\frac{t}{T_0}$

In addition, we introduce a normalized amplitude U as:

$A(z,\tau)=\sqrt{P_0}e^{\frac{- \alpha Z}{2}}U(z,\tau)$

We now take into consideration a space scale normalized to the fiber length as:

$\zeta=\frac{z}{L}$

Where L is the fiber length. Thus, it turns out that $U(\zeta,\tau)$ satisfies

$\frac{\partial U}{\partial \zeta}=+i \frac{L}{2 L_D} sgn(\beta_2) \frac{\partial^2 U}{\partial \tau^2}-i e^{- \alpha Z \zeta} \frac{L}{L_{NL}} (|U|^2 U)$

Where $sgn(\beta_2)=\pm 1$ depending on the sign of the coefficient $\beta_2$ and

$L_D=\frac{T_0 ^ 2}{|\beta_2|}$

$L_{NL}=\frac{1}{\gamma_0 P_0}$

The concept behind this normalization process is to exclude any kind of overflow error that may occur during solving the PDE with finite difference method. The derivation of the dimensionless transformation for the nonlinear Schrödinger equation can be found [here](https://arxiv.org/pdf/1612.00358). The plotting functions originally came from [here](https://github.com/OleKrarup123/NLSE-vector-solver). The 3-Stage Explicit Fractional Order Runge-Kutta Method algorithm in Mathematica can be found [here](https://www.mdpi.com/2504-3110/7/3/245).