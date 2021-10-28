import math

try:
    from gmpy2 import root
except ImportError:
    def root(x, n):
        return x ** (1 / n)

try:
    from math import prod
except ImportError:
    def prod(xs):
        return math.exp(sum(map(math.log, xs)))
