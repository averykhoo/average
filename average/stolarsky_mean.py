import math


def stolarsky_mean_1(xs, r, s):
    """
    https://www.mdpi.com/2227-7390/8/8/1320
    page 2

    note that the special cases of the 2-var stolarsky mean are no longer the same special cases for this
    """
    n = len(xs)
    assert n != 0
    assert r != 0
    assert s != 0
    assert s != r
    prod_xs = math.prod(xs)
    return ((r ** 2 / s ** 2) *
            (sum(x ** (n * s) for x in xs) - n * (prod_xs ** s)) /
            (sum(x ** (n * r) for x in xs) - n * (prod_xs ** r))
            ) ** (1 / (n * (s - r)))


def stolarsky_mean_2(ws, xs, ys, r, s):
    """
    https://www.mdpi.com/2227-7390/8/8/1320
    page 3
    renamed A to W since those seem to be weights, and `as` is a python reserved word)
    no clue what Y is supposed to be

    note that the special cases of the 2-var stolarsky mean are no longer the same special cases for this

    note that if you set W = 1 and Y = 0, then you can sort of get the lehmer mean / gini mean
    """
    assert len(ws) == len(xs) == len(ys)
    assert r != 0
    assert s != 0
    assert s != r

    return ((r ** 2 / s ** 2) *
            sum(w * (x ** s - y ** s) ** 2 for w, x, y in zip(ws, xs, ys)) /
            sum(w * (x ** r - y ** r) ** 2 for w, x, y in zip(ws, xs, ys))
            ) ** (1 / (2 * (s - r)))

if __name__ == '__main__':

    xs = [2, 3, 4, 5, 6, 7, 8, 9]
    for i in range(10):
        print(10*i, stolarsky_mean_1([10*i+x for x in xs],1,2))