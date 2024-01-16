# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: data_ex2.py
"""
import numpy as np 
import math

# Domain definition:
x1 = 0
y1 = 0
x2 = 1
y2 = 1
h = 0.1

# Plot options:
plot_solution = 'y'
plot_error = 'n'

def f(x, y):
    return 0

def g1(x,y):
    return math.exp(x+y)

def g2(x,y):
    return math.exp(x+y)

def g3(x,y):
    return math.exp(x+y)

def g4(x,y):
    return math.exp(x+y)

def mu(x, y):
    return 1

def beta(x, y):
    b0 = 0
    b1 = 0
    return [b0, b1]

def sigma(x,y):
    return 2

def u_an(x,y):
    return math.exp(x+y)