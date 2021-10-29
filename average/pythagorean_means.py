import math
import operator
from functools import reduce


def arithmetic_mean(xs):
    return sum(xs) / len(xs)


def harmonic_mean(xs):
    return len(xs) / sum(1 / x for x in xs)


# `math.prod` is only available in Python 3.8 and above
if hasattr(math, 'prod'):
    def geometric_mean(xs):
        return math.prod(xs) ** (1 / len(xs))
else:
    def geometric_mean(xs):
        try:
            return math.exp(arithmetic_mean(list(map(math.log, xs))))
        except ValueError:
            return reduce(operator.mul, xs, 1) ** (1 / len(xs))
