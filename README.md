# Modelling microalgae growth

The code from a research project I worked on which aimed to develop a mathematical model for simulating microalgae growth for CO<sub>2</sub> utilisation. 

Some modelling was attempted using <a href="http://www.daetools.com">DAE Tools</a>, which is an "equation-based, object-oriented modelling, simulation and optimisation" software based in Python. However, the main optimisation results were generated using the sequential least-squares quadratic programming (SLSQP) algorithm within the SciPy optimize method.

The model is first solved for an initial set of parameters and the sum of squared errors (SSE) between the model solution and some experimental data is found. The SSE is then minimised to find the optimised parameter values which best fit the model to the data.
