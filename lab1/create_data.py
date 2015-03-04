#!/usr/bin/python
#coding=UTF-8

__author__ = 'fedoramy'

from sympy import *
import numpy as np
import copy
import re
from Scientific.Geometry import Vector


class Create_Data():
    def input_from_file(self, path="input.txt"):
        lines = []
        try:
            f = open(path, 'r')
            lines = f.readlines()
        finally:
            f.close()
        return lines

    def str_to_float(self, str):
        try:
            return float(str)
        except:
            return 0

    def str_to_int(self, str):
        try:
            return int(str)
        except:
            return 0

    def string_to_float(self, str):
        try:
            x = []
            for i in range(len(str)):
                x.append(float(str[i]))
            return x
        except:
            return 0

    def input_data(self):
        lines = self.input_from_file()
        arguments = lines[1].split()
        x_string = lines[2].split()
        x0 = self.string_to_float(x_string)
        eps1 = float(lines[3])
        eps2 = float(lines[4])
        eps3 = float(lines[5])
        iters = int(lines[6])
        return [lines[0], arguments, x0, eps1, eps2, eps3, iters]

    def is_empty_list(self, l):
        if len(l) == 0:
            return true
        else:
            return false

    def create_matrix_A(self, function, arguments):
        matrix = []
        mas = []
        n = len(arguments)
        # print(n)
        for i in range(n):
            for j in range(n):
                if i == j:
                    p = re.compile('([+-]?\d+|[+-]?\d+\.\d+)\*?' + arguments[i] + '\*{2}2')
                    try:
                        mas.append(float(p.findall(function)[0]))
                    except:
                        mas.append(0.0)
                else:
                    p = re.compile('([+-]?\d+|[+-]?\d+\.\d+)\*?' + arguments[i] + '\*' + arguments[j])
                    z = re.compile('([+-]?\d+|[+-]?\d+\.\d+)\*?' + arguments[j] + '\*' + arguments[i])
                    p_l = p.findall(function)
                    if self.is_empty_list(p_l):
                        p_l.append(0.0)

                    z_l = z.findall(function)
                    if self.is_empty_list(z_l):
                        z_l.append(0.0)
                    try:
                        mas.append((float(p_l[0]) + float(z_l[0]))/2)
                    except:
                        mas.append(0.0)

            matrix.append(copy.deepcopy(mas))
            del mas[:]
        Matrix = np.copy(matrix)
        Matrix = 2*Matrix
        return Matrix

    def create_b(self, function, arguments):
        n = len(arguments)
        b = []
        for i in range(n):
            p = re.compile('([+-]?\d+|[+-]?\d+\.\d+)\*?' + arguments[i]+'[^*]')
            p_l = p.findall(function)
            if self.is_empty_list(p_l):
                b.append([0.0])
            else:
                b.append([float(p_l[0])])
        B = np.copy(b)
        return B

class Alg_Quadratic(object):

    MAX_ITER =100000

    class ToFile(object):
        def __init__(self, i, xk, f, f_deriv, norm1, norm2, norm3):
            self.i = i
            self.xk = xk
            self.f = f
            self.f_deriv = f_deriv
            self.norm1 = norm1
            self.norm2 = norm2
            self.norm3 = norm3

        def tostr(self):
            rep ="######################################################"+'iters: '+ str(self.i)+'\nXi='\
                 +str(self.xk)+'\nF(Xi)='+str(self.f)+"\nF'(Xk)="\
                  + str(self.f_deriv)+'\n||Xi+1 - X||='+ str(self.norm1)+'\n||F(Xi+1)-F(Xi)||='\
                + str(self.norm2)+"\n||F'(Xi+1)-F'(Xi)||="+str(self.norm3) + '\n'
            return rep

    def __init__(self, matrix, vector, x0, e1, e2, e3, iterations):
        self.A = matrix
        self.B = vector
        self.X0 = x0
        self.eps1 = e1
        self.eps2 = e2
        self.eps3 = e3
        self.iters = iterations

    def function(self, x):
        '''
         return 1/2(Ax, x) + (b, x)
        '''
        return 1.0/2 * float(np.dot(np.transpose(np.dot(self.A, x)), x)) + float(np.dot(np.transpose(self.B), x))

    def dif(self, x):
        return np.dot(self.A, x) + self.B

    def Hk(self, x):
        return -1*(self.dif(x))

    def alpha(self, x):
        '''
        :param x: input porint
        :return: -(f'(x), h)/(Ah,h)
        '''
        m = np.dot(np.transpose(self.dif(x)), self.Hk(x))
        n = np.dot(np.transpose(np.dot(self.A, self.Hk(x))) , self.Hk(x))
        return -1*float(m)/float(n)

    def norm_max(self, x):
        return np.linalg.norm(x, np.inf)

    def iteration(self, path = 'output.txt'):
        xk = self.X0
        filename = open(path, 'a+')
        i = 0
        while(true):
            i += 1
            if i>self.MAX_ITER:
                print('iter more than '+str(self.MAX_ITER))
                filename.close()
                return 0
            xk1 = xk + self.alpha(xk)*self.Hk(xk)
            f = self.function(xk)
            f_deriv = self.dif(xk)
            f1 = self.function(xk1)
            f_deriv1 = self.dif(xk1)
            norm1 = self.norm_max(xk1 - xk)
            norm2 = abs(f1 - f)
            norm3 = self.norm_max(f_deriv1)
            a = self.ToFile(i, xk1, f1, f_deriv1, norm1, norm2, norm3)
            self.output(filename, a)
            if norm1 < self.eps1 and (norm2 < self.eps2) and (norm3 < self.eps3):
                filename.close()
                print(a.tostr())
                return 0
            xk = xk1

    def output(self, f, a):
        f.write(a.tostr())



ob = Create_Data()
list1 = ob.input_data()
matrix = ob.create_matrix_A(list1[0], list1[1])
b = ob.create_b(list1[0], list1[1])
X = np.matrix(list1[2])
x = np.transpose(X)
eps1 = list1[3]
eps2 = list1[4]
eps3 = list1[5]
iters = list1[6]

print(matrix)
print(b)
print(x)

obj = Alg_Quadratic(matrix, b, x, eps1, eps2, eps3, iters)
obj.iteration()











