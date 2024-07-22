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
h = 0.01

# Problem definition:
E = 10E+9
A = 1
P = 100

# Plot options:
plot_solution = 'y'
plot_error = 'y'

def f(x):
    return 0

def g1(x):
    return 0

def g2(x):
    return 1

def h1(x):
    return 0

def h2(x):
    return P

def mu(x):
    return E*A

def beta(x):
    return 0

def sigma(x):
    return 0

def u_analytical(x):
    return (P/E*A)*x

def grad_u_analytical(x):
    return P/E*A