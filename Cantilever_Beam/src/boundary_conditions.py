# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE FEM 1D
Title: Boundary_condition.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

def impose_boundary_conditions(A, F, mesh):
    g = np.zeros(2*mesh.ndof)
    A, F, g = impose_Dirichlet("left", A, F, g, mesh)
    #A, F, g = impose_Dirichlet("right", A, F, g, mesh)
    F = impose_Neumann("right", F, mesh)
    return A, F, g


def impose_Dirichlet(wall, A_, F_, g_, mesh_):
    bc_index_list = []
    g = g_
    g_wall = np.zeros(2*mesh_.ndof)

    if wall == "left":
        bc_index_list.extend([0, 1])
        g_wall[0] = data.g1(mesh_.nodes[0])
        g_wall[1] = data.g2(mesh_.nodes[0])
        g[:2] = g_wall[:2] 

    elif wall == "right":
        bc_index_list.extend([-1, -2])
        g_wall[-1] = data.g3(mesh_.nodes[-1])
        g_wall[-2] = data.g4(mesh_.nodes[-1])
        g[-2:] = g_wall[-2:] 

    # Lifting:
    F = F_ - A_.dot(g_wall)

    A = A_
    for bc_index in bc_index_list:
        for i in range(len(A_)):
            A[bc_index][i] = 0
            A[i][bc_index] = 0
        A[bc_index][bc_index] = 1
        F[bc_index] = 0 
    
    return A, F, g

def impose_Neumann(wall, F_, mesh_):
    F = F_

    if wall == "left":
        F[0] = F_[0] - data.h1(mesh_.nodes[0])
        F[1] = F_[1] - data.h2(mesh_.nodes[0])

    elif wall == "right":
        F[-2] = F_[-2] + data.h3(mesh_.nodes[-1])
        F[-1] = F_[-1] + data.h4(mesh_.nodes[-1])

    return F