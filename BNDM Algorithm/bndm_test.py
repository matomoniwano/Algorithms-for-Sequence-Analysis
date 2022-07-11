
from bndm import build_nfa_and, bndm
import numpy as np

def test_bndm():
    sequence = b"ACGGGCTAGCTACGACGTACGATCAGCT"
    P = b"AGCT"
    nfa = build_nfa_and(P[::-1])
    NRESULTS = 5
    results = np.zeros(NRESULTS, dtype=np.uint64)
    nresults = bndm(*nfa, sequence, results)
    assert nresults == 2
    assert results[0] == 10
    assert results[1] == 27
    assert results[2] == 0
    assert results[3] == 0
    assert results[4] == 0
