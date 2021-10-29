import math
import statistics


def arithmetic_geometric_mean(xs):
    """
    https://en.wikipedia.org/wiki/arithmetic-geometric_mean
    """
    a_mean = statistics.mean(xs)
    g_mean = statistics.geometric_mean(xs)

    while abs(g_mean - a_mean) / max(g_mean, a_mean) > 1e-15:
        a_next = (a_mean + g_mean) / 2
        g_next = math.sqrt(g_mean * a_mean)
        a_mean = a_next
        g_mean = g_next

    return g_mean


def arithmetic_harmonic_mean(xs):
    """
    https://mathworld.wolfram.com/Arithmetic-HarmonicMean.html
    """
    return statistics.geometric_mean(xs)


def geometric_harmonic_mean(xs):
    """
    https://en.wikipedia.org/wiki/geometric-harmonic_mean
    """
    g_mean = statistics.geometric_mean(xs)
    h_mean = statistics.harmonic_mean(xs)

    while abs(g_mean - h_mean) / max(g_mean, h_mean) > 1e-15:
        g_next = math.sqrt(g_mean * h_mean)
        h_next = 2 / (1 / g_mean + 1 / h_mean)
        g_mean = g_next
        h_mean = h_next

    return g_mean
