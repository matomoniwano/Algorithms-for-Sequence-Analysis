import argparse

import numpy as np

from naive_pattern_matching import naive_pattern_matching, get_patterns, get_text

def test_naive_pattern_matching():
    T = "TTACGTATTTTTCGAGTACGTT"
    Ps = ["ACGT", "TTTT", "ACGATAT"]
    found_sol = [True, True, False]
    n_sol = [2,2,0]
    pos_sol = [[2,17], [7,8], []]
    comp_sol = [27, 36, 21]

    for i, P in enumerate(Ps):
        found, n, pos, comp = naive_pattern_matching(P,T)
        assert found == found_sol[i], f"Pattern: {P}"
        assert n == n_sol[i], f"Pattern: {P}"
        assert pos == pos_sol[i], f"Pattern: {P}"
        assert comp == comp_sol[i], f"Pattern: {P}"


def test_get_pattern():
    p = argparse.ArgumentParser(description="DNA naive pattern matching")
    pat = p.add_mutually_exclusive_group()
    pat.add_argument("-P", "--pattern",
        help="immediate pattern to be matched")
    pat.add_argument("-p", "--patternfile",
        help="name of file containing patterns (one per line)")

    args = p.parse_args()
    args.pattern = "ACCGGATA"
    assert len(get_patterns(args)) == 1

    args.pattern = None
    args.patternfile = "assignment01_patterns.txt"
    Ps = get_patterns(args)
    print(type(Ps[0]))
    if isinstance(Ps[0], bytes):
        for i in range(len(Ps)):
            Ps[i] = Ps[i].decode('ASCII')
    assert not isinstance(Ps[0], bytes)
    print(Ps)
    assert len(Ps) == 5
    assert Ps[0] == "ACGATTAGCTA"
    assert Ps[1] == "ACC"
    assert Ps[2] == "ATGTCAT"
    assert Ps[3] == ""
    assert Ps[4] == "ACTAGATC"

def test_get_text():
    p = argparse.ArgumentParser(description="DNA naive pattern matching")
    txt = p.add_mutually_exclusive_group()
    txt.add_argument("-T", "--text",
        help="immerdiate text to be searched")
    txt.add_argument("-t", "--textfile",
        help="name of file containing text (will be read in one piece)")

    args = p.parse_args()
    t = "CAGTCAGTCATCATGCGTATCAGCTGATCTATCGGGCGCGCGCGTATCATC"
    args.text = t

    text = get_text(args)
    if isinstance(text, bytes): text=text.decode('ASCII')
    assert text == t

    args.text = None
    args.textfile = "assignment01_text.txt"
    text = get_text(args)
    if isinstance(text, bytes): text=text.decode('ASCII')
    assert text == t
