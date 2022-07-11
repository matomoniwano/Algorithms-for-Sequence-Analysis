import argparse

import numpy as np

def naive_pattern_matching(P, T):
    """
    Implementation of the naive patern matching algorithm.
    
    Input:
    P: Pattern (string)
    T: Text in which the pattern is to be searched (string)
    Output:
    found: boolean if the pattern was found or not
    n: Number of occurrences of the pattern
    positions: List of all positions a pattern starts
    comps: Number of all comparisons made
    """

    found = False
    n = 0
    positions = []
    comps = 0

   # variables for length of input string
    p_len = len(P) 
    t_len = len(T)

    #if the length of pattern is longer than the length of text then exit
    if p_len > t_len:
      return 1
    
    # for each character
    for i in range(t_len - p_len + 1):
      if T[i:i+p_len] == P: # if the pattern match then
        found = True # make it true
        n = n + 1 # and count
        positions.append(i) # add matched position
        comps = comps + p_len # comparison for full match is always the length of pattern
      else: # if the first character is not matched
        counter = 0 # this is a counter for comparison
        while(True): # while loop breaks when there is a mismatch from the comparison
          if T[i + counter] == P[counter]:
            comps = comps + 1
            counter = counter + 1
          else:
            comps = comps + 1
            break

    return found, n, positions, comps

def get_text(args):
    """
    Input:
    args: argument parser
    Output: Text (string)
    """

    text = ""
    # if input argument is a string
    if args.text and not args.textfile:
      text = args.text
      text.strip()
    # if input arugment is a text file
    elif args.textfile and not args.text:
      file = open(args.textfile)
      for i in file.readlines():
        text+=str(i.strip())
 
    
    return text


def get_patterns(args):
    """
    Input:
    args: argument parser
    Output: List of all patterns (list of strings)
    """

    Ps = []
    # if input argument is a string
    if args.pattern:
      Ps.append(args.pattern.strip()) 
    # if input argument is a text file
    elif args.patternfile:
      file = open(args.patternfile)
      for i in file.readlines():
        Ps.append(i.strip()) # append in the list


    return Ps


def main(args):
    T = get_text(args)
    Ps = get_patterns(args)

    for P in Ps:  # iterate over patterns
        if len(P) == 0: continue  # skip empty patterns
        print(f"> {P}")
        found, n , positions, comps = naive_pattern_matching(P, T)
        print(f"Pattern {P}:\n Found: {found}\n Occurred: {n} times \n Positions: {positions}\n Comparisons: {comps}")


def get_argument_parser():
    p = argparse.ArgumentParser(description="DNA naive pattern matching")
    pat = p.add_mutually_exclusive_group(required=True)
    pat.add_argument("-P", "--pattern",
        help="immediate pattern to be matched")
    pat.add_argument("-p", "--patternfile",
        help="name of file containing patterns (one per line)")
    txt = p.add_mutually_exclusive_group(required=True)
    txt.add_argument("-T", "--text",
        help="immerdiate text to be searched")
    txt.add_argument("-t", "--textfile",
        help="name of file containing text")
    return p


if __name__ == "__main__":
    main(get_argument_parser().parse_args())
