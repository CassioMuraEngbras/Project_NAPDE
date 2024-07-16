# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: data_ex2.py

Example: Bar with axial loads. 
"""
import numpy as np 
import math

#https://www.wolframalpha.com/input?i2d=true&i=-+D%5Bu%2C%7Bx%2C2%7D%5D+%2B+u+%3D+sin%5C%2840%29x%5C%2841%29%5C%2844%29+u%5C%2840%290%5C%2841%29%3D0%5C%2844%29+u%5C%2840%292+pi%5C%2841%29+%3D+0

# Domain definition:
x1 = 0
x2 = 2
h = 0.01

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
    return 100

# Material properties definition:
E = 10E+9
A = 1
def mu(x):
    return E*A

def beta(x):
    return 0

def sigma(x):
    return 0

def u_analytical(x):
    return 0.5*math.sin(x)

def grad_u_analytical(x):
    return 0.5*math.cos(x)