import math
import statistics

from average.pythagorean_means import geometric_mean


def generalized_mean(xs, dim=1):
    """
    Also known as the "power mean"
    https://en.wikipedia.org/wiki/Generalized_mean

    special cases:
    dim=-inf -> minimum
    dim=-1   -> harmonic mean
             -- geometric-harmonic mean would fit in here
    dim=0    -> geometric mean
             -- logarithmic mean would fit in here (arithmetic-logarithmic-geometric mean inequality)
             -- arithmetic-geometric mean would fit in here
    dim=1    -> arithmetic mean
    dim=2    -> root mean square (quadratic mean)
    dim=3    -> cubic mean
             -- contra-harmonic mean would fit in here
    dim=inf  -> maximum
    """
    # extreme special cases
    if dim == -math.inf:
        return min(xs)
    if dim == math.inf:
        return max(xs)

    # special cases
    if dim == -1:
        return statistics.harmonic_mean(xs)
    if dim == 0:
        try:
            return statistics.geometric_mean(xs)
        # occasionally the statistics way of calculating exp(mean(map(log, xs))) doesn't quite work as intended
        # eg. because you can't take the log of zero
        except statistics.StatisticsError:
            return geometric_mean(xs)
    if dim == 1:
        return statistics.mean(xs)
    if dim == 2:
        return root_mean_square(xs)

    # technically only defined over positive numbers
    if dim < 0:
        if any(x <= 0 for x in xs):
            raise statistics.StatisticsError('values must be positive when dim is negative')
    elif dim % 2:
        if any(x < 0 for x in xs):
            raise statistics.StatisticsError('negative values not allowed when dim is not even')

    # general case
    return statistics.mean(x ** dim for x in xs) ** (1 / dim)


def root_mean_square(xs):
    """
    special case since (excepting some really smart optimizations) I think x * x should be faster than x ** 2
    and I assume math.sqrt is somehow better than x ** 0.5
    """
    return math.sqrt(statistics.mean(x * x for x in xs))


def generalized_f_mean(xs, f, f_inv=None):
    """
    https://en.wikipedia.org/wiki/Quasi-arithmetic_mean

    The function `f` must be monotonic, and `f_inv` must be its inverse.
    It's not a great idea to compute the inverse of some arbitrary function f, but it can be done, albeit slowly
    Newton's method might work better here, but binary search will suffice for now

    returns f_inv(mean([f(x) for x in xs]))
    """
    if f_inv is None:
        def f_inv(y):
            hi = 1
            lo = -1

            # find bounds
            if f(1) < f(2):  # monotonically increasing
                while f(hi) < y:
                    hi *= 2
                while f(lo) > y:
                    lo *= 2
            else:
                while f(hi) > y:
                    hi *= 2
                while f(lo) < y:
                    lo *= 2

            # binary search
            while lo <= hi and hi - lo > 1e-15:
                x = (lo + hi) / 2
                if f(x) < y < f(hi) or f(x) > y > f(hi):
                    lo = x
                elif f(x) < y < f(lo) or f(x) > y > f(lo):
                    hi = x
                else:
                    return x
            return lo

    return f_inv(statistics.mean(map(f, xs)))


def _generalized_mean(xs, dim):
    """
    illustrates that the generalized f-mean is a further generalization of the generalized mean
    """
    if dim == 0:
        return _geometric_mean(xs)
    else:
        return generalized_f_mean(xs, lambda x: math.pow(x, dim), lambda y: math.pow(y, 1 / dim))


def _geometric_mean(xs):
    """
    illustrates how much of a special case the geometric mean is, since it doesn't really work with (dim := 0)

    effectively computes the following:
    ```python
    product = 1
    for value in data:
        product *= value

    geometric_mean = product ** (1 / len(data))
    ```
    """
    return generalized_f_mean(xs, math.log, math.exp)


def _log_sum_exp(xs):
    """
    https://en.wikipedia.org/wiki/LogSumExp
    """
    return generalized_f_mean(xs, math.log, math.exp) + math.log(len(xs))


def _p_norm(xs, p):
    """
    https://en.wikipedia.org/wiki/Lp_space#The_p-norm_in_finite_dimensions
    """
    if p <= 0:
        raise ValueError('p must be positive')
    if p == math.inf:
        return max(map(abs, xs))
    return generalized_f_mean(map(abs, xs), lambda x: x ** p, lambda y: y ** (1 / p))
