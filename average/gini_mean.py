def gini_mean(xs, r, s):
    """
    https://www.jstor.org/stable/2322316

    notice that if you set s = 1, then you get the r-th root of the lehmer mean
    the gini mean mentioned in https://www.mdpi.com/2227-7390/8/8/1320 is equivalent if you squint the right way

    :param xs:
    :param r:
    :param s:
    :return:
    """
    return (sum(x ** (r + s) for x in xs) / sum(x ** r for x in xs)) ** (1 / r)
