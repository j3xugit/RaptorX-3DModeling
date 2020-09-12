import numpy as np
import math
import sys

def LogNormal(x, a, b):
    x0 = np.array(x) + sys.float_info.epsilon
    B = 2*b*b
    y = np.exp(-(np.log(x0)-a)**2/B)
    y = y /math.sqrt(2*math.pi) / b
    return y

def OneLogNormal(x, a, b, c, d):
    y = d - c *LogNormal(x, a, b)
    return y

def TwoLogNormal(x, a0, b0, c0, a1, b1, c1, d):
    y = d - c0* LogNormal(x, a0, b0) -c1 * LogNormal(x, a1, b1)
    return y

def ThreeLogNormal(x, a0, b0, c0, a1, b1, c1, a2, b2, c2, d):
    y = d - c0* LogNormal(x, a0, b0) -c1 * LogNormal(x, a1, b1) - c2 * LogNormal(x, a2, b2)
    return y
