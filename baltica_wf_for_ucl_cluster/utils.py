from itertools import combinations


def all_against_all_design(conditions):
    return [(a, b) for a, b in combinations(conditions, 2)]
