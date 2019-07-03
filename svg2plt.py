import numpy as np
from matplotlib import pyplot as plt
import re
import xml.etree.ElementTree as ET
import sys

N = 10
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

def linear_bezier(t, pts):
    return (1 - t) * pts[0] + t * pts[1]

def quadratic_bezier(t, pts):
    return (1 - t) ** 2 * pts[0] + \
        2 * t * (1 - t) * pts[1] + \
        t ** 2 * pts[2]

def cubic_bezier(t, pts):
    return (1 - t) ** 3 * pts[0] + \
        3 * t * (1 - t) ** 2 * pts[1] + \
        3 * t ** 2 * (1 - t) * pts[2] + \
        t ** 3 * pts[3]

bezier = {'l': linear_bezier,
          'q': quadratic_bezier,
          'c': cubic_bezier}

orders = {'l': 1, 'q': 2, 'c': 3}

t = np.linspace(0, 1, N)

initx, inity = None, None
lastx, lasty = 0, 0

path = parse_path(path_string)
Xs = np.array([])
Ys = np.array([])

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
                Xs = np.append(Xs, func(t, xs[i:i+order+1]))
                Ys = np.append(Ys, func(t, ys[i:i+order+1]))
                i += order
        else:
            raise Exception('Unsupported cmd ' + cmd)
        if initx is None:
            initx, inity = lastx, lasty
    elif len(p) == 1:
        Xs = np.append(Xs, linear_bezier(t, [lastx, initx]))
        Ys = np.append(Ys, linear_bezier(t, [lasty, inity]))
        lastx, lasty = initx, inity
        initx, inity = None, None

plt.plot(Xs, Ys, 'k.')
plt.axis('equal')
plt.show()
