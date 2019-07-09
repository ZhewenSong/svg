import numpy as np
from matplotlib import pyplot as plt
import re
import xml.etree.ElementTree as ET
import sys
from sympy import symbols, integrate, exp, I, pi, sin, cos

t = symbols("t")

supported_cmd = 'mlqcz'

inputfile = sys.argv[1]
path_string = ''

tree = ET.parse(inputfile)
for ele in tree.iter():
    if ele.tag.endswith('path'):
        path_string += ele.get('d')

def parse_path(path_string):
    path = []
    coords = ''
    for c in path_string:
        if c.isalpha():
            if path and coords:
                data = list(map(float, coords.split()))
                if data:
                    path[-1].append(data)
            path.append([c])
            coords = ''
        else:
            coords += c
    if path and coords:
        data = list(map(float, coords.split()))
        if data:
            path[-1].append(data)
    return path

def linear_bezier(pts, dt=0):
    return (1 - (t - dt)) * pts[0] + (t - dt) * pts[1]

def quadratic_bezier(pts, dt=0):
    return (1 - (t - dt)) ** 2 * pts[0] + \
        2 * (t - dt) * (1 - (t - dt)) * pts[1] + \
        (t - dt) ** 2 * pts[2]

def cubic_bezier(pts, dt=0):
    return (1 - (t - dt)) ** 3 * pts[0] + \
        3 * (t - dt) * (1 - (t - dt)) ** 2 * pts[1] + \
        3 * (t - dt) ** 2 * (1 - (t - dt)) * pts[2] + \
        (t - dt) ** 3 * pts[3]

bezier = {'l': linear_bezier,
          'q': quadratic_bezier,
          'c': cubic_bezier}

orders = {'l': 1, 'q': 2, 'c': 3}

initx, inity = None, None
lastx, lasty = 0, 0

path = parse_path(path_string)
Xs = []
Ys = []

dt = 0
for p in path:
    if len(p) == 2:
        cmd, coords = p
        if cmd == 'M':
            lastx, lasty = coords[0], coords[1]
        elif cmd == 'm':
            lastx += coords[0]
            lasty += coords[1]
        elif cmd.lower() in 'lqc':
            xs = [lastx] + coords[::2]
            ys = [lasty] + coords[1::2]
            order = orders[cmd.lower()]
            if cmd.islower():
                for i in range(1, len(xs)):
                    xs[i] += lastx
                    ys[i] += lasty
                    if i % order == 0:
                        lastx, lasty = xs[i], ys[i]
            i = 0
            while i + order < len(xs):
                func = bezier[cmd.lower()]
                Xs.append((func(xs[i:i+order+1], dt), (t, dt, dt+1)))
                Ys.append((func(ys[i:i+order+1], dt), (t, dt, dt+1)))
                i += order
                dt += 1
        else:
            raise Exception('Unsupported cmd ' + cmd)
        if initx is None:
            initx, inity = lastx, lasty
    elif len(p) == 1:
        if lastx != initx or lasty != inity:
            Xs.append((linear_bezier([lastx, initx], dt), (t, dt, dt+1)))
            Ys.append((linear_bezier([lasty, inity], dt), (t, dt, dt+1)))
            dt += 1
            lastx, lasty = initx, inity
        initx, inity = None, None


# for X, Y in zip(Xs, Ys):
#     print("{{If[{}<=t<={},{},0], If[{}<=t<={},{},0]}},"
#         .format(X[1][1], X[1][2], str(X[0]).replace('**', '^'), 
#             Y[1][1], Y[1][2], str(Y[0]).replace('**', '^')))
# N = 10
# coeff = []
# for k in range(-N, N+1):
#     integrals = []
#     for X, Y in zip(Xs, Ys):
#         integrals.append("NIntegrate[({} + I*{}) * Exp[-2Pi*I*{}*t/{}], {{{}, {}, {}}}]"
#             .format(str(X[0]).replace('**', '^'), str(Y[0]).replace('**', '^'), 
#         k, dt, *X[1]))
    
#         coeff.append('({}) * Exp[2Pi*I*{}*t/{}]'.format('+'.join(integrals), k, dt))
# print('+'.join(coeff))

# N = 10
# an = []
# bn = []
# freq = range(-N, N+1)
# for k in freq:
#     an.append(sum(integrate(X * cos(2*pi*I*k*t/dt), limit) for X, limit in Xs) / dt)
# print(coef)
# plt.plot(Xs, Ys, 'k.')
# plt.axis('equal')
# plt.show()
