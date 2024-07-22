# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: data_ex2.py

Example: Bar with axial loads. 
"""
import numpy as np 
import math

# Domain definition:
x1 = 0
x2 = 2
h = 0.1

# Problem definition:
L = x2 - x1
E = 2*np.power(10, 6)
I = 1.666/np.power(10, 5)
P = 2000

# Plot options:
plot_solution = 'y'
plot_error = 'n'

def f(x):
    return 0

## Add the variables y1, theta1, y2, theta2.
def g1(x):
    return 0

def g2(x):
    return 0

#def g3(x):
#    return 0

#def g3(x):
#    return 0

#def h1(x):
#    return 0

#def h2(x):
#    return 0

def h3(x):
    return P

def h4(x):
    return 0

def mu(x):
    return E*I

def beta(x):
    return 0

def sigma(x):
    return 0

def u_analytical(x):
    return -P/(6*E*I)*np.power(x, 3) + P*L/(2*E*I)*np.power(x, 2)

def grad_u_analytical(x):
    return -P/(2*E*I)*x*x + P*L/(E*I)*x

def M_analytical(x):
    return P*(x - L)

def V_analytical(x):
    return P + 0*x