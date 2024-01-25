from algebra import LCSgraph
from algebra.variants import parse_hgvs
from algebra.extractor.local_supremal import local_supremal
from bisect import bisect, bisect_left


def interferes(start, end, flat_supremal):
    i = bisect(flat_supremal, start)
    j = bisect_left(flat_supremal, end)
    if i == j and i % 2 == 0:
        return False
    else:
        return True

def interference_type(start, end, flat_supremal):
    bisect_a = bisect(flat_supremal, start)
    bisect_b = bisect_left(flat_supremal, end)
    if bisect_a == bisect_b:
        if bisect_a % 2 != 0:
            return "within local supremal"
    elif bisect_a == bisect_b - 1:
        return "one boundary crossed"
    else:
        return "spanning one or multiple local supremal"


def rna(reference, variant, intervals):
    graph = LCSgraph.from_variant(reference, parse_hgvs(variant))
    # print("\n".join(to_dot(reference, graph, labels=False)))
    print(graph.supremal)
    local_supremals = local_supremal(reference, graph)
    print("-------\nlocal supremals:")
    for local in local_supremals:
        print(local.start, local.end)
    print("".join([" " + b for b in reference]))
    flat_supremals = [a for b in [(l.start, l.end) for l in local_supremals] for a in b]
    print(flat_supremals)
    for (a, b) in intervals:
        if a == b:
            raise ValueError("Interval start equal to interval end.")
        if interferes(a, b, flat_supremals):
            print(a, b, interference_type(a, b, flat_supremals))
        else:
            print(a, b, "OK")


def main():
    reference = "CTCTAGAGACTTTATTTTCCACG"
    variant = "[1_6delinsGTCTC;14A>C;17_18insA;21A>C]"
    intervals = [(8, 9), (9, 13), (17, 18), (8, 10), (8, 17)]

    rna(reference, variant, intervals)


if __name__ == "__main__":
    main()