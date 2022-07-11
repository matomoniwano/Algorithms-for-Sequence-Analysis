import argparse  # for command line interface

def build_nfa_and(P):
    """Build an NFA from a pattern"""
    """
    Input:
        Pattern P (String)
    Output:
        mask: Dictionary of all masks for all different characters in the pattern (dict)
        accept:
    """

    masks = dict()
    accept_state = 0

    """
    TODO: For a given pattern P calculate all bit masks and store them in the dictionary masks.
    Also return the accepting state.
    
    """
    
    # Get all the unique characters from given pattern
    unique = list(set(P))
    unique.sort()
    
    # bit length is determined based on the length of given pattern string
    m = len(P)
    accept_bit = m * "0"
      
    #Go over the pattern and create masks for each character    
    for i in unique:
      mask = ""
      for j in P:
        if i == j:
          mask = mask + "1" # if matched, then change bit to 1
        else:
          mask = mask + "0"
      mask = mask[::-1]
      masks[i] = int(mask, 2)
    
    #else mask is created when pattern is searched in a text that contains character that is not present in pattern (just in case)
    else_mask = m * "0"
    masks["else"] = int(else_mask, 2) # create a key called "else" which is mapped to 0
    accept_bit = list(accept_bit)
    # accept state is 1 followed by 0
    accept_bit[0] = "1"
    accept_state_bit = ''.join(accept_bit)
    accept_state = int(accept_state_bit,2)
    
    return masks, accept_state



def build_nfa_or(P):
    """
    Input:
        Pattern P (String)
    Output:
        mask: Dictionary of all masks for all different characters in the pattern (dict)
        accept:
    """

    masks = dict()
    accept_state = 0

    """
    TODO: For a given pattern P calculate all bit masks and store them in the dictionary masks.
    Also return the accepting state.
    Note: To build the or nfa you need to invert the bit logic
    """
    unique = list(set(P))
    unique.sort()
    m = len(P)
    accept_bit = m * "1" # accept state for "OR" is 0 followed by 1
      
    # repeat the same algorithm from AND implementation but in reverse
    for i in unique:
      mask = ""
      for j in P:
        if i == j:
          mask = mask + "0"
        else:
          mask = mask + "1"
      mask = mask[::-1]
      masks[i] = int(mask, 2)
    
    else_mask = m * "1"
    masks["else"] = int(else_mask, 2)
    accept_bit = list(accept_bit)
    accept_bit[0] = "0"
    accept_state_bit = ''.join(accept_bit)
    accept_state = int(accept_state_bit,2)
    
    return masks, accept_state


def shift_and(masks, accept, text, N):
    """
    Input:
        masks: Dictionary of all masks (dict)
        accept: Which state is the accepting state (int)
        text: Text (String)
        N: Maximum number of positions to store in results
    Output:
        results: List of all positions at which a pattern ends (list)
        k: number of matches (int)
    """

    k = 0
    results = []

    """
    TODO:
    Implement the Shift-And algorithm as descirbed in the lecture.
    Masks and accept are computed by the build_nfa_and function and provided as parameters.
    The results list should only include the first N results.
    k is the number of all occurrences.
    """
    len_t = len(text)
    len_p = len(format(accept, "b")) # get the length of given pattern length
    status = len_p * "0" # status is used to keep track of current state
    get_bin = lambda x, n: format(x, 'b').zfill(n) # a lambda function to convert int to binary
    goal = format(accept, "b")
    
    for i in range(len(text)):
      match = False
      status_list = list(status)
      status_list[len_p - 1] = "1" #00001
      if(text[i] in masks):
        i_th_mask = list(get_bin(masks[text[i]], len_p)) #mask
      else:
        i_th_mask = list(get_bin(masks["else"], len_p))
      for bit in range(len(i_th_mask)): #shift algorithm
        if status_list[0] == "1" and status_list[0] == i_th_mask[0]:
          match = True
        if  status_list[bit] == "1" and i_th_mask[bit] == "0":
          status_list[bit] = "0"
        elif status_list[bit] == "1" and i_th_mask[bit] == "1":
            status_list[bit - 1] = "1"
            status_list[bit] ="0"
        
            
        
      status = ''.join(status_list)
      
      #we have a new status and now we will check for accept state
      if match != False:
        k = k + 1
        if len(results) < N:
          results.append(i)
    return k, results


def shift_or(mask, accept, text, N):
    """
    Input:
        masks: Dictionary of all masks (dict)
        accept: Which state is the accepting state (int)
        text: Text (String)
        N: Maximum number of positions to store in results
    Output:
        results: List of all positions at which a pattern ends (list)
        k: number of matches (int)
    """

    k = 0
    results = []

    """
    TODO:
    Implement the Shift-Or algorithm.
    This one is similiar to the shift and algorithm but you need to invert the bit logic.
    Masks and accept are computed by the build_nfa_or function and provided as parameters.
    The results list should only include the first N results.
    k is the number of all occurrences.
    
    """
    get_bin = lambda x, n: format(x, 'b').zfill(n) # a lambda function to convert int to binary
    len_t = len(text)
    len_p = len(format(mask["else"], "b")) # get the length of given pattern length
    status = len_p * "1" # status is used to keep track of current state
    goal = format(accept, "b")
    
    for i in range(len(text)):
      match = False
      status_list = list(status)
      status_list[len_p - 1] = "0" #00001
      if(text[i] in mask):
        i_th_mask = list(get_bin(mask[text[i]], len_p)) #mask
      else:
        i_th_mask = list(get_bin(mask["else"], len_p))
      for bit in range(len(i_th_mask)):
        if status_list[0] == "0" and status_list[0] == i_th_mask[0]:
          match = True
        if  status_list[bit] == "0" and i_th_mask[bit] == "1":
          status_list[bit] = "1"
        elif status_list[bit] == "0" and i_th_mask[bit] == "0":
            status_list[bit - 1] = "0"
            status_list[bit] ="1"
        
            
        
      status = ''.join(status_list)
      
      #we have a new status and now we will check for accept state
      if match != False:
        k = k + 1
        if len(results) < N:
          results.append(i)
    
    return k, results



def get_text(args):
    if args.text is not None:
        return args.text
    with open(args.textfile, "r") as ftext:
        text = ftext.read()
    return text


def get_patterns(args):
    if args.pattern is not None:
        return [args.pattern]  # list with single item
    with open(args.patternfile, "r") as fpat:
        Ps = [pattern.strip() for pattern in fpat.readlines()]
    return Ps


def main(args):
    alg = args.algorithm
    T = get_text(args)  # bytes object
    Ps = get_patterns(args)  # list of bytes objects
    build_nfa = build_nfa_and if alg == "and" else build_nfa_or
    find_matches = shift_and if alg == "and" else shift_or
    NRESULTS = args.maxresults
    for P in Ps:  # iterate over patterns
        if len(P) == 0: continue  # skip empty patterns
        nfa = build_nfa(P)
        nresults, results = find_matches(*nfa, T, NRESULTS)
        if nresults > NRESULTS:
            print("! Too many results, showing first {NRESULTS}")
            nresults = NRESULTS
        print(*list(results[:nresults]), sep="\n")


def get_argument_parser():
    p = argparse.ArgumentParser(description="DNA Motif Searcher")
    pat = p.add_mutually_exclusive_group(required=True)
    pat.add_argument("-P", "--pattern",
        help="immediate pattern to be matched")
    pat.add_argument("-p", "--patternfile",
        help="name of file containing patterns (one per line)")
    txt = p.add_mutually_exclusive_group(required=True)
    txt.add_argument("-T", "--text",
        help="immerdiate text to be searched")
    txt.add_argument("-t", "--textfile",
        help="name of file containing text (will be read in one piece)")
    p.add_argument("-a", "--algorithm", metavar="ALGORITHM",
        default="and", choices=("and", "or"),
        help="algorithm to use ('and' (default), 'or')")
    p.add_argument("--maxresults", "-R", type=int, default=10_000,
        help="maximum number of results to show (10_000)")
    return p


if __name__ == "__main__":
    main(get_argument_parser().parse_args())
