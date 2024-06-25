# Welcome to ODE Docs

This is the homework 4 from Computational Physics course, done by Daniel Lacayo Zúñiga (B84180). Is important to mention that this project is inspired on `calculator`project used on the article [Build Your Python Project Documentation With MkDocs](https://realpython.com/python-project-documentation-with-mkdocs/), according to the instructions.

This site contains the project documentation for the `ode` project.
Its aim is to give us a framework to build our project documentation using Python, MkDocs, mkdocstrings, and the Material for MkDocs theme.

In this homework, we work with three numerical methods well-known on computational fields; Euler's method, Runge-Kutta 2 and Runge-Kutta 4.

## Euler

The Euler's method can give us a good aproximation depending of the problem and the amout of iterations that we need in our solution. Generally, with this method, a calculus that is twice more precise require the twice of computational resources.

Relevant equation:

$$
result[i] = result[i - 1] + h * function(x_0 + (i - 1) * h, result[i - 1])
$$

## RK2

The Runge-Kutta method is in fact a family of different order methods which bring us a better aproximation without the necessity to consider higher orders in the Taylor's expansion of Euler's method. We want to avoid this last point, because of the fact that is complicated to know the derivate of the function which we are evaluating in the right side of the ODE.

The idea of RK2 method is to use the middle point in order to evaluate the Euler's method. Euler's method is applied in the point $$t$$ to evaluate the derivate to approximate the function in the point $$x = t + h;$$ for the other hand RK2 method use the middle point $$t + h/2,$$ obtaining a better aproximation for the same value of $$h.$$

Relevant equation:

$$
result[i] = result[i - 1] + 0.5 * (k_1 + k_2)
$$

## RK4

The last metodology can applied even to more points between $$x(t)$$ and $$x(t + h)$$, making Taylor's expansions. In that case, it can to group terms of order $$h^3, h^4, $$ etc.; to cancell those expresions.

The problem to do this is that expresions involve more complicated when we increase the aproximation order. In general terms, the rule is that $$4^{\\rm t_o}$$ order correspond to the best compromise between complexity and aproximation error. This method is the most used commonly in order to resolve ODEs.

Relevant equation:

$$
result[i] = result[i - 1] + (k_1 + 2 * kh_2 + 2 * k_3 + k_4) / 6
$$

## Project layout


    mkdocs.yml        # The configuration file.
    docs/
        index.md      # Introduction to three numerical methods, with relevant equations.
        reference.md  # Documentation of four functions.


## Project Overview

::: ode
