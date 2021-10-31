import math

from functools import lru_cache


@lru_cache(maxsize=65536)
def log_mean(*xs):
    """
    generalized logarithmic mean
    https://en.wikipedia.org/wiki/Logarithmic_mean#Generalization

    in base 2, log_mean(1, 2) == 1, which breaks the math
    not sure if there's a similar problem in base e, but doesn't seem like it
    """
    xs = sorted(xs)  # helps with caching and dealing with duplicates
    if len(xs) == 0:
        raise TypeError('expected at least 1 arguments, got 0')
    elif xs[0] < 0:
        raise ValueError('cannot log a negative number')
    elif xs[0] == 0:
        # as xs[0] -> 0, log_mean -> 0, even though log(0) is undefined
        return 0
    elif xs[0] == xs[-1]:
        # by definition, min <= mean <= max, so if min == max, then min == mean
        # this also takes care of single-argument means
        return xs[0]
    elif len(xs) == 2:
        # https://en.wikipedia.org/wiki/Logarithmic_mean
        ret = (xs[1] - xs[0]) / (math.log(xs[1]) - math.log(xs[0]))
        assert xs[0] < ret < xs[-1], (xs, ret)
        return ret
    else:
        # https://www.survo.fi/papers/logmean.pdf (page 7, formula 14)
        # todo: does this match the output of the generalized stolarsky mean?
        ret = (len(xs) - 1) * (log_mean(*xs[1:]) - log_mean(*xs[:-1])) / (math.log(xs[-1]) - math.log(xs[0]))
        assert xs[0] < ret < xs[-1], (xs, ret)
        return ret
