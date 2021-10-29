import statistics


def lehmer_mean(xs, dim):
    """
    https://en.wikipedia.org/wiki/Lehmer_mean

    special cases:
    dim=-inf -> minimum
    dim=0    -> harmonic mean
    dim=0.5  -> geometric mean
    dim=1    -> arithmetic mean
    dim=2    -> contra-harmonic mean
    dim=inf  -> maximum
    """
    # extreme special cases
    if dim == -float('inf'):
        return min(xs)
    if dim == float('inf'):
        return max(xs)

    # special cases
    if dim == 0:
        return statistics.harmonic_mean(xs)
    if dim == 0.5:
        return statistics.geometric_mean(xs)
    if dim == 1:
        return statistics.mean(xs)
    if dim == 2:
        return contra_harmonic_mean(xs)

    #
    if any(x <= 0 for x in xs):
        raise statistics.StatisticsError('negative values not allowed when -inf < dim <= 0')
    return sum(x ** dim for x in xs) / sum(x ** (dim - 1) for x in xs)


def weighted_lehmer_mean(xs, dim, weights=None):
    if weights is None:
        return lehmer_mean(xs, dim)
    return sum(w * (x ** dim) for w, x in zip(weights, xs)) / sum(w * (x ** (dim - 1)) for w, x in zip(weights, xs))


def contra_harmonic_mean(xs):
    return sum(x * x for x in xs) / sum(xs)
