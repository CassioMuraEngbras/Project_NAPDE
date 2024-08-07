# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: local_element.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data
import src.quadrature as quadrature

class Local_element:
    def __init__(self, point_1, point_2, point_3, nodes_indexes):
        self.point_1 = point_1
        self.point_2 = point_2
        self.point_3 = point_3
        self.nodes_indexes = nodes_indexes

        x1, y1 = point_1[0], point_1[1]
        x2, y2 = point_2[0], point_2[1]
        x3, y3 = point_3[0], point_3[1]
        self.J = [x2 - x1, y2 - y1],[x3 - x1, y3 - y1]
        self.T = [[1, x1, y1],[1, x2, y2], [1, x3, y3]]

    def phi(self, n, x_, y_):
        # Coefficients of the basis function:
        b = [[1, 0, 0], [0, 1, 0], [0, 0, 1]][n]
        a0, a1, a2 = np.linalg.solve(self.T , b)
        return a0 + a1*x_ + a2*y_

    def dphi(self, n, x_, y_):
        # Coefficients of the basis function:
        b = [[1, 0, 0], [0, 1, 0], [0, 0, 1]][n]
        a0, a1, a2 = np.linalg.solve(self.T , b)
        return [a1, a2]

    def A_local(self):
        A_local_ = np.zeros([3, 3])
        K_local = np.zeros([3, 3]) # Local stiffness matrix
        V_local = np.zeros([3, 3]) # Local advection matrix
        M_local = np.zeros([3, 3]) # Local mass matrix
        det_J = np.linalg.det(self.J)
        for m in range(3):
            for n in range(3):
                for quadrature_weight, quadrature_point in zip(quadrature.weights(), quadrature.points()):
                    x_, y_ = quadrature.map(self, quadrature_point[0], quadrature_point[1])
                    phi_m, phi_n = self.phi(m, x_, y_), self.phi(n, x_, y_)
                    dphi_m, dphi_n = self.dphi(m, x_, y_), self.dphi(n, x_, y_)

                    K_local[m][n] += abs(det_J)*quadrature_weight*(data.mu(x_, y_)*np.dot(dphi_m, dphi_n))
                    V_local[m][n] += abs(det_J)*quadrature_weight*(np.dot(data.beta(x_, y_), dphi_n)*phi_m)
                    M_local[m][n] += abs(det_J)*quadrature_weight*(data.sigma(x_, y_)*phi_m*phi_n)
                    
        A_local_ = K_local + V_local + M_local
        return A_local_

    def F_local(self):
        F_local_ = np.zeros(3)
        det_J = np.linalg.det(self.J)
        for m in range(3):
            for quadrature_weight, quadrature_point in zip(quadrature.weights(), quadrature.points()):
                x_, y_ = quadrature.map(self, quadrature_point[0], quadrature_point[1])
                phi_m = self.phi(m, x_, y_)

                F_local_[m] += abs(det_J)*quadrature_weight*(data.f(x_, y_)*phi_m)
        return F_local_