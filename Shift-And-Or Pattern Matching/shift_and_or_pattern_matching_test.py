from shift_and_or_pattern_matching import build_nfa_and, build_nfa_or, shift_and, shift_or

def test_build_nfa_and():
    masks, accept = build_nfa_and("ABABC")
    assert masks["A"] == int("00101", 2)
    assert masks["B"] == int("01010", 2)
    assert masks["C"] == int("10000", 2)

def test_build_nfa_or():
    masks, accept = build_nfa_or("ABABC")
    assert masks["A"] == int("11010", 2)
    assert masks["B"] == int("10101", 2)
    assert masks["C"] == int("01111", 2)

def test_shift_and():
    masks, accept = build_nfa_and("ABABC")
    text = "ABABABCCACABAC"
    k, results = shift_and(masks, accept, text, 1)
    assert k == 1
    assert results[0] == 6


    masks, accept = build_nfa_and("ABAB")
    text = "ABABABABABABABABABABABABABABABAB"
    k, results = shift_and(masks, accept, text, 5)
    assert k == 15
    assert len(results) == 5
    assert results == [3,5,7,9,11]

def test_shift_or():
    masks, accept = build_nfa_or("ABABC")
    text = "ABABABCCACABAC"
    k, results = shift_or(masks, accept, text, 1)
    assert k == 1
    assert results[0] == 6

    masks, accept = build_nfa_or("ABAB")
    text = "ABABABABABABABABABABABABABABABAB"
    k, results = shift_or(masks, accept, text, 5)
    assert results == [3,5,7,9,11]
    assert k == 15
    assert len(results) == 5
