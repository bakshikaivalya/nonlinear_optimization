#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A primal-dual convex hull algorithm to obtain the non-redundant constraints of a given n dimensional linear constraint
inequality Ax < b. Nonlinear interior point solvers commonly used in planning control, may have a very high number of
linear constraints which may grow linearly with the time complexity. Unloading all constraints into the solver leads to
a high dimensional Hessian which grows quadratically with the number of decision variables as well as constraints
(which become slack decision variables in interior point methods). The algorithm and 2d example is an implementation of
a fast and inexpensive convex hull approach to identifying the non-redundant constraints i.e. redundant constraint
removal.

Author: Kaivalya Bakshi
@Date: 15 Mar 2019
"""

import numpy as np
import numpy.matlib
from scipy import linalg
from scipy.spatial import ConvexHull
import scipy.optimize
import matplotlib.pyplot as plt
import random
from scipy.optimize import minimize

class reduceConstr():
    def __init__(self, A, b):  # _TODO_
        self.A = A
        self.b = b

    def linConstrInteriorScore(self, point):  # score = max(A*x - b) is required to be negative for x to be in the interior of the feasible region
        A = self.A
        b = self.b
        score = np.max(np.append((A @ point.reshape(A.shape[1], 1) - b).reshape(b.size, 1), [0]))  #np.max([0, np.max(A @ point - b)])  # linear constraint Ax - b < 0 # np.max(np.max(A @ point - b))  #
        return score

    def getInteriorPoint(self):
        A = self.A
        b = self.b
        c = linalg.lstsq(A, b)[0]  # least squares solution to Ax = b expected to be atleast in the vicinity of the interior of the feasible region
        if (np.max(list((A @ c - b).reshape(b.size, ))) > 0):
            copt = minimize(self.linConstrInteriorScore, x0=c, method='nelder-mead',
                            options={'xtol': 1e-8, 'disp': True}).x.reshape(A.shape[1], 1)  # 'inner-most interior point using Nelder-Mead method
            print('Least squares solution is not an interior point of the feasible region.')
            print('Interior score max(0, vec(A@copt-b)) is ' + str(np.max(A @ copt - b)) + ' at point ' + str(copt))
        else:
            copt = c
            print('Least squares solution is an interior point of the feasible region.')
            if (np.max(list((A @ copt - b).reshape(b.size, ))) > 0):
                print(
                    'Could not find an interior point of the feasible region using Nelder-Mead method. Try something else.')
                raise SystemExit
        return copt

    def nonredundantCstr(self, copt):
        A = self.A
        b = self.b
        b_translated = b - (A @ copt).reshape(b.shape[0], 1)  # scaling to make copt the origin
        D = A / b_translated  # element wise division to obtain the normalized form D*x + ones < 0. Note that points in the matrix D are the dual vertices corresponding to primal facets/constraint planes
        try:
            hull = ConvexHull(D)  # dual of the primal feasible region polytope obtained as the convex hull of dual points of constraints
        except:
            print('Could not obtain dual polytope of primal feasible region using the convex hull approach.')
            raise SystemExit
        hull_ind = np.unique(hull.simplices.flat)
        DD_transpose = D @ D.transpose()
        A_red = np.zeros((len(hull_ind), 2))
        b_red = np.zeros((len(hull_ind), 1))
        for red_ind in range(len(hull_ind)):
            A_red[red_ind][:] = A[red_ind][:]
            b_red[red_ind][:] = b[red_ind][:]
        return D, hull_ind, hull, A_red, b_red

    def x1x2vectorsLine2d(self, xvecinput, a1, a2, b):  # function for plotting 2d lines by returning vector of x2 or x1 values given vector of values in an interval of x1 or x2
        if (np.abs(a2) < 10**-4) and (np.abs(a1) > 10**-4):  # if a2 \approx 0
            x1vec = b * np.ones(xvecinput.shape) / a1
            x2vec = xvecinput
        elif (np.abs(a2) > 10**-4) and (np.abs(a1) < 10**-4):  # if a1 \approx 0
            x1vec = xvecinput
            x2vec = b * np.ones(xvecinput.shape) / a2
        else:  # both a1, a2 are not \approx 0
            x1vec = xvecinput
            x2vec = (b * np.ones(xvecinput.shape) - a1 * xvecinput) / a2
        return x1vec, x2vec

if __name__ == '__main__':
    # 2d example of obtaining non-redundant constraints given the constraint matrices in Ax < b
    A = np.array([[1, 2], [4, 2], [-1, 1], [0, -1], [-1, 0], [1, 1], [-2, 1]])
    b = np.array([[4], [12], [1], [0], [0], [5], [5]])
    constrReductionExample = reduceConstr(A, b)
    copt = constrReductionExample.getInteriorPoint()
    D, hull_ind, hull, A_red, b_red = constrReductionExample.nonredundantCstr(copt)
    # plotting the 2d example
    x1vec = np.linspace(-1, 4, 100)
    plt.plot(copt[0], copt[1], 'k+', markersize=5, label='interior point (primal space)')
    plt.plot(D[:, 0], D[:, 1], 'ko', markerSize=10, label='dual vertices (dual space)')
    plt.plot(D[hull.vertices, 0], D[hull.vertices, 1], 'ro--', label='convex hull (dual space)')
    for i in range(A.shape[0]):
        x1, x2 = constrReductionExample.x1x2vectorsLine2d(x1vec, A[i][0], A[i][1], b[i])
        if (i == min(hull_ind)):
            label = 'original constraints (primal space)'
        else:
            label = '_no_label_'
        plt.plot(x1, x2, 'b', label=label, linewidth=2)
        if i in hull_ind:
            x1, x2 = constrReductionExample.x1x2vectorsLine2d(x1vec, A[i][0], A[i][1], b[i])
            if (i==min(hull_ind)):
                label = 'non-redundant constraints (primal space)'
            else:
                label = '_no_label_'
            plt.plot(x1, x2, 'r--', linewidth=1, label=label)
    plt.grid()
    plt.title(str(len(hull_ind)) + ' non-redundant constraints (primal space)')
    plt.legend()
    plt.xlabel('x1')
    plt.ylabel('x2')