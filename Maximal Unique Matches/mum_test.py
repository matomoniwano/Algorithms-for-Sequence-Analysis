from MUM import count_mums, compute_pos_manber_myers, compute_lcp
def test_mum():
    T = b"miississippii"
    S = b"mississippi"
    n1, n2 = len(T), len(S)  # lengths of the individual genomes
    T = bytes(T+S)
    n = len(T)
    pos = compute_pos_manber_myers(T)
    lcp = compute_lcp(T, pos)
    nmums, length = count_mums(T, pos, lcp, n1, minlen=0, show=False)
    assert nmums == 2
    assert length == 12
