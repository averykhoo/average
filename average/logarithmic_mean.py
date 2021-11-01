import math

from functools import lru_cache


@lru_cache(maxsize=65536)
def log_mean(*xs):
    """
    generalized logarithmic mean
    https://en.wikipedia.org/wiki/Logarithmic_mean#Generalization

    there are issues with numerical stability
    it likely does NOT match the output of the generalized stolarsky mean

    in base 2, log_mean(1, 2) == 1, which breaks the math
    not sure if there's a similar problem in base e, but doesn't seem like it

    ERROR when xs = [92, 93, 94, 95, 96, 97, 98, 99]
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
        ret = (xs[1] - xs[0]) / (math.log(xs[1] / xs[0]))
        assert xs[0] < ret < xs[-1], (xs, ret)
        return ret
    else:
        # https://www.survo.fi/papers/logmean.pdf (page 7, formula 14)
        ret = (len(xs) - 1) * (log_mean(*xs[1:]) - log_mean(*xs[:-1])) / math.log(xs[-1] / xs[0])
        assert xs[0] < ret < xs[-1], (xs, ret)
        return ret


def _numerically_unstable_logarithmic_mean(xs):
    """
    https://www.survo.fi/papers/logmean.pdf (page 3, formula 4)
    """

    @lru_cache(maxsize=0xFFFF)
    def fib(n):
        """
        fibonacci sequence
        """
        if n == 0:
            return 0
        elif n == 1:
            return 1
        else:
            return fib(n - 1) + fib(n - 2)

    elems = []
    for i, x in enumerate(xs):
        elems.append(x / math.prod((math.log(x) - math.log(y)) for j, y in enumerate(xs) if i != j))
    print(elems)
    return fib(len(xs) - 1) * sum(elems)


def stolarsky_mean(xs, r, s):
    """
    extending means of two variables to several variables
    page 10, formula 8.3

    :param xs:
    :param r:
    :param s:
    :return:
    """
    return (log_mean([x ** s for x in xs]) / log_mean([x ** r for x in xs])) ** (1 / (s - r))


if __name__ == '__main__':
    # print(log_mean(*range(1,100)))
    # print(log_mean(92, 93, 94, 95, 96, 97, 98, 99))
    # print(log_mean(52, 53, 54, 55, 56, 57, 58, 59))
    # print(log_mean(2, 3, 4, 5, 6, 7, 8, 9))
    xs = [2, 3, 4, 5, 6, 7, 8, 9]
    for i in range(10):
        print(10 * i, log_mean(*[10 * i + x for x in xs]))
    # print(log_mean( 73, 74, 75, 76, 77, 78,79))
