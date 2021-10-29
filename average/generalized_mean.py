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
    return math.sqrt(statistics.mean(x * x for x in xs))
