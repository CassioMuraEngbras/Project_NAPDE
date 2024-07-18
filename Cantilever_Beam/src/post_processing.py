# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE FEM 1D
Title: post_processing.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

''' Basis function:
def psi(mesh, x):
    h = mesh.h
    psi1 = 1 - 3*np.power(x/h, 2) + 2*np.power(x/h, 3)
    psi2 = -x*np.power(1 - x/h, 2)
    psi3 = np.power(x/h, 2)*(3 - 2*x/h)
    psi4 = h*np.power(x/h, 2)*(1 - x/h)
    return np.array([psi1, psi2, psi3, psi4])

def dpsi(mesh, x):
    h = mesh.h
    dpsi1 = -(6/h)*(x/h)*(1 - x/h)
    dpsi2 = -1 + 4*x/h - 3*np.power(x/h, 2)
    dpsi3 = -dpsi1
    dpsi4 = (x/h)*(2 - 3*x/h)
    return np.array([dpsi1, dpsi2, dpsi3, dpsi4])
'''
def ddpsi(mesh, x):
    h = mesh.h
    ddpsi1 = -(6/np.power(h, 2))*(1 - 2*x/h)
    #ddpsi2 = -(2/h)*(3*x/h - 2)
    ddpsi2 = (2/h)*(3*x/h - 2)
    ddpsi3 = - ddpsi1
    #ddpsi4 = (2/h)*(1 - 3*x/h)
    ddpsi4 = -(2/h)*(1 - 3*x/h)
    return np.array([ddpsi1, ddpsi2, ddpsi3, ddpsi4])

def dddpsi(mesh, x):
    h = mesh.h
    dddpsi1 = 12/np.power(h, 3)
    #dddpsi2 = -6/np.power(h, 2)
    dddpsi2 = 6/np.power(h, 2)
    dddpsi3 = -12/np.power(h, 3)
    #dddpsi4 = -6/np.power(h, 2)
    dddpsi4 = 6/np.power(h, 2)
    return np.array([dddpsi1, dddpsi2, dddpsi3, dddpsi4])

def bending_moment(mesh, U):
    M = []
    for i in range(len(mesh.coord) - 1):
        U_element = np.array([U[2*i], U[2*i+1], U[2*i+2], U[2*i+3]])
        M.append(-data.E*data.I*np.dot(ddpsi(mesh, 0), U_element))
        if i == len(mesh.coord) - 2:
           M.append(-data.E*data.I*np.dot(ddpsi(mesh, mesh.h), U_element)) 
    return M

def shear_force(mesh, U):
    V = []
    for i in range(len(mesh.coord) - 1):
        U_element = np.array([U[2*i], U[2*i+1], U[2*i+2], U[2*i+3]])
        V.append(-data.E*data.I*np.dot(dddpsi(mesh, 0), U_element))
        if i == len(mesh.coord) - 2:
           V.append(-data.E*data.I*np.dot(dddpsi(mesh, mesh.h), U_element)) 
    return V
    
'''' -------------- Error Calculation -----------------
# Should have been called interpolate:
def u_h(x , x1, y1, x2, y2):
        return (x - x1)*(y2 - y1)/(x2 - x1) + y1

def err(x, x1, y1, x2, y2):
    return (data.u_analytical(x) - u_h(x, x1, y1, x2, y2))**2

def err_grad(x, x1, y1, x2, y2):
    return (data.grad_u_analytical(x) - u_h(x, x1, y1, x2, y2))**2

def compute_gradient(U, mesh):
    gradient_U = np.zeros(mesh.ndof)
    for i in range(len(U) - 1):
        gradient_U[i] = (U[i+1] - U[i])/mesh.h
    # Extrapolation for computing the last element:
    gradient_U[-1] = 2*gradient_U[-2] - gradient_U[-3] 

    return gradient_U

def compute_errors(U, mesh):
    # Gaussian quadrature points:
    quad_points = [-math.sqrt(3/5), 0 , + math.sqrt(3/5)]
    quad_weights = [5/9, 8/9 , 5/9]

    # Computation of the gradient:
    grad_U = compute_gradient(U, mesh)

    U_analytical = np.zeros(mesh.ndof)
    grad_U_analytical = np.zeros(mesh.ndof)
    for i in range(mesh.ndof):
        U_analytical[i] = data.u_analytical(mesh.nodes[i])
        grad_U_analytical[i] = data.grad_u_analytical(mesh.nodes[i])

    L2_error = 0
    for element in mesh.elements:
        for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
            x_map = 0.5*(element.x2 - element.x1)*quadrature_point + 0.5*(element.x2 + element.x1)
            L2_error += 0.5*(element.x2 - element.x1)*quadrature_weight*err(x_map, element.x1, U[element.index_node_1], element.x2, U[element.index_node_2])
    L2_error = math.sqrt(L2_error)

    H1_semi_error = 0
    for element in mesh.elements:
        for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
            x_map = 0.5*(element.x2 - element.x1)*quadrature_point + 0.5*(element.x2 + element.x1)
            H1_semi_error += 0.5*(element.x2 - element.x1)*quadrature_weight*err_grad(x_map, element.x1, grad_U[element.index_node_1], element.x2, grad_U[element.index_node_2])
    H1_semi_error = math.sqrt(H1_semi_error)

    H1_error = math.sqrt(L2_error**2 + H1_semi_error**2)

    if data.plot_error == 'y':
        fig = plt.figure()
        plt.scatter(mesh.coord, U, color ='black', marker = 'o', label = "Finite Element Solution")
        plt.scatter(mesh.coord, U_analytical, color ='purple', marker = 'x', label = "Analytical Solution")
        plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
        plt.title('Finite Element Solver')
        plt.xlabel('x')
        plt.ylabel('u(x)')
        plt.figtext(0.15, 0.83, "L2-error = " + '%.6f' %L2_error)
        plt.show()
        
    return L2_error, H1_error
    '''