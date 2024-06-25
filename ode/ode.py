# ode/ode.py

"""Provide several sample numeric math methods calculations.

This module allows the user to make mathematical calculations.

Examples:
    >>> from ode import ode.py
    >>> ode.function(2, 4)
    2

The module contains the following functions:

- `function(x, y)` - Returns a number.
- `euler(x0, xf, fprima, ite)` - Apply the Euler's numerical method.
- `rk2(x0, xf, fprima, ite)` - Apply RK2 numerical method.
- `rk4(x0, xf, fprima, ite)` - Apply RK4 numerical method.
"""

import numpy as np

def function(x, y):
    """
    A function that bring us a number.
    
    Examples:
        >>> function(2, 3)
        2
        >>> function(5, 0)
        5

    Args:
        x (float): Independent variable.
        y (float): Dependent variable.

    Returns:
        float: The value computed by the function at (x, y).
    """
    return x

def euler(x0, xf, function, ite):
    """
    Solves an ODE using the Euler method.

    Examples:
    	>>> euler(-2, 2, function, 40)
    	[-2, -2.2, -2.39, -2.57, -2.74, -2.9, -3.05, -3.19, -3.32, -3.44, -3.55, -3.65, -3.74, -3.82, -3.89, -3.95
    	-4, -4.04, -4.07, -4.09, -4.1, -4.1, -4.09, -4.07, -4.04, -4, -3.95, -3.89, -3.82, -3.74, -3.65, -3.55, -3.44
    	-3.32, -3.19, -3.05, -2.9, -2.74, -2.57, -2.39]

    Args:
        x0 (float): Initial value of x.
        xf (float): Final value of x.
        function (function): Function of two variables.
        ite (int): Number of iterations.

    Returns:
        np.ndarray: Array with the solution of the ODE at the given points.
    """
    h = (xf - x0) / ite
    result = np.zeros(ite)
    result[0] = x0
    for i in range(1, ite):
        result[i] = result[i - 1] + h * function(x0 + (i - 1) * h, result[i - 1])
    return result

def rk2(x0, xf, function, ite):
    """
    Solves an ODE using the second-order Runge-Kutta method (RK2).

    Examples:
        >>> rk2(-2, 2, function, 40)
        [-2, -2.1975, -2.385, -2.5625, -2.73, -2.8875, -3.035, -3.1725, -3.3, -3.4175, -3.525, -3.6225, -3.71, -3.7875,
        -3.855, -3.9125, -3.96, -3.9975, -4.025, -4.0425 -4.05, -4.0475, -4.035, -4.0125, -3.98, -3.9375, -3.885,
        -3.8225, -3.75, -3.6675, -3.575, -3.4725, -3.36, -3.2375, -3.105, -2.9625, -2.81, -2.6475, -2.475, -2.2925]

    Args:
        x0 (float): Initial value of x.
        xf (float): Final value of x.
        function (function): Function of two variables.
        ite (int): Number of iterations.

    Returns:
        np.ndarray: Array with the solution of the ODE at the given points.
    """
    h = (xf - x0) / ite
    result = np.zeros(ite)
    result[0] = x0
    for i in range(1, ite):
        k1 = h * function(x0 + (i - 1) * h, result[i - 1])
        k2 = h * function(x0 + (i - 0.5) * h, result[i - 1] + 0.5 * k1)
        result[i] = result[i - 1] + 0.5 * (k1 + k2)
    return result

def rk4(x0, xf, function, ite):
    """
    Solves an ODE using the fourth-order Runge-Kutta method (RK4).
    
    Examples:
        >>> rk4(0, 1, function, 10)
        [-2, -2.195, -2.38, -2.555, -2.72, -2.875, -3.02, -3.155, -3.28, -3.395, -3.5, -3.595, -3.68, -3.755, -3.82,
        -3.875, -3.92, -3.955, -3.98, -3.995, -4, -3.995, -3.98, -3.955, -3.92, -3.875, -3.82, -3.755, -3.68, -3.595,
        -3.5, -3.395, -3.28, -3.155, -3.02, -2.875, -2.72, -2.555, -2.38, -2.195]

    Args:
        x0 (float): Initial value of x.
        xf (float): Final value of x.
        function (function): Function of two variables.
        ite (int): Number of iterations.

    Returns:
        np.ndarray: Array with the solution of the ODE at the given points.
    """
    h = (xf - x0) / ite
    result = np.zeros(ite)
    result[0] = x0
    for i in range(1, ite):
        k1 = h * function(x0 + (i - 1) * h, result[i - 1])
        k2 = h * function(x0 + (i - 0.5) * h, result[i - 1] + 0.5 * k1)
        k3 = h * function(x0 + (i - 0.5) * h, result[i - 1] + 0.5 * k2)
        k4 = h * function(x0 + i * h, result[i - 1] + k3)
        result[i] = result[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return result
