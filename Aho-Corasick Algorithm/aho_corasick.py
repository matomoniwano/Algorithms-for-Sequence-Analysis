import collections   # has a deque for breadth-first-search
import argparse

#####################################################
# Aho-Corasick

class ACNode():
    
    def __init__(self, parent=None, letter=None, depth=0, label=""):
        self.targets = dict()  # children of this node
        self.lps = None        # lps link of this node
        self.parent = parent   # parent of this node
        self.letter = letter   # letter between parent and this node
        self.out = []          # output function of this node
        if parent is None:
            self.depth = depth # number of chars from root to this n ode
            self.label = label # string from root to this node
        else:
            self.depth = parent.depth + 1
            self.label = parent.label + letter

        
    def delta(self, a):
        """
        Transit to next node upon processing character a.
        INPUT:
        - self is the object. you can acces parameters like parent using self.parent
        - a: text character
        OUTPUT:
        - q: the next state (another ACNode instance)
        TODO: Implement the delta function for the current node:
        What is the target of self after reading character a?
        """
        
        if a in self.targets:
          q = self.targets[a]
        else:
          q = ACNode()

        return q


    def bfs(self):
        """yields each node below and including self in BFS order"""
        Q = collections.deque([self])
        while len(Q) > 0:
            node = Q.popleft()
            yield node
            Q.extend(node.targets.values())

    def __str__(self):
        s = ""
        for (i,node) in enumerate(self.bfs()):
            lps = node.lps.label if node.lps!=None else "[none]"
            s += "{0}: label={1}, targets={2}, lps={3}, output={4}\n".format(
                i, node.label, list(node.targets.keys()), lps, node.out)
        return s


def AC_build(P):
    """build AC autmaton for list of patterns P, return its root node."""
    # Build a root for the trie
    root = ACNode()
    # Build the trie, pattern by pattern
    for (i,p) in enumerate(P):
        node = root
        for a in p:
            if a in node.targets:
                node = node.targets[a]
            else:
                newnode = ACNode(parent=node, letter=a)
                node.targets[a] = newnode
                node = newnode
        node.out.append(i)
    # Walk through trie in BFS-order to build lps
    for node in root.bfs():
        if node.parent is None: continue
        node.lps = node.parent.lps.delta(node.letter) if node.depth>1 else root
        node.out.extend(node.lps.out)
    return root


def search_with_AC(P, T):
    """
    INPUT:
    - List of Patterns P
    - Text T
    OUTPUT:
    - yield each triple (start, stop, pattern_index)
    """

    root = AC_build(P)

    # TODO: Implement the pattern search using Aho-Corasick-Algorithm.
    # The Trie is built in AC_build.
    for i in range(len(T)): # looping through each character from given text
      # if the character is leading to the next node
      # then move the node to that node
      if T[i] in root.targets:
        root = root.delta(T[i]) 
      # if the character is not leading to the next node then there are several options
      else:       
        #this option is when it is not leading back to the root, instead finding the longest matching suffix
        if root.lps is not None: 
            if root.lps.label != "":
              root = root.lps
            else:
              label = root.label
              next = True
              while next != False:
                temp = AC_build(P)
                label = label[1:]
                for z in label:
                  temp = temp.delta(z)
                if temp.lps is not None:
                  root = temp
                  break
                if len(label) == 0:
                  root = AC_build(P)
                  break
        # going back to root and reset.
        elif root.lps is None: 
          label = root.label
          next = True
          while next != False:
            temp = AC_build(P)
            label = label[1:]
            for l in label:
              temp = temp.delta(l)
            if temp.lps is not None:
              root = temp
              break
            if len(label) == 0:
              root = AC_build(P)
              break 
        if T[i] in root.targets: #comparing with the current letter
            root = root.delta(T[i])
            
      #Checking if output has any value then yield
      if len(root.out) != 0:
        for output in root.out:
          p_len = len(P[output])
          start = i - p_len + 1
          end = i + 1
          pattern_index = output
          found = (start, end, pattern_index)
          print(found)
          if (T[start:end] == P[output]):
            yield found

    #  When text is finished iterating (checking to see the remaining matching patter from last terminal state)
    done = True
    label_last = root.label
    while done != False:
        temp_last = AC_build(P)
        label_last = label_last[1:]
        if len(label_last) == 0:
          break
        else:
          for letter in label_last:
              temp_last = temp_last.delta(letter)
          if temp_last.lps is not None:
              root = temp_last      
              break
            
    if len(root.out) != 0:
        for output in root.out:
            p_len = len(P[output])
            start = i - p_len + 1
            end = i + 1
            pattern_index = output
            found = (start, end, pattern_index)
            if (T[start:end] == P[output]):
                yield found


def main(args):
    T = args.text
    P = args.pattern
    ret = search_with_AC(P,T)
    print(list(ret))

def get_argument_parser():
    p = argparse.ArgumentParser(description="DNA Motif Searcher")
    p.add_argument("-P", "--pattern", required=True, nargs="+",
        help="immediate pattern to be matched")
    p.add_argument("-T", "--text", required=True,
        help="immerdiate text to be searched")
    return p

if __name__ == "__main__":
    main(get_argument_parser().parse_args())
