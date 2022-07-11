# Suffixarrays und Suffixbaeume
###########################################################################
import argparse

###########################################################################
# SUFFIXTREE

#Node class
class node:
    
    #initialize the class
    def __init__(self):
        self.remainder = 0 #remainder counter
        self.activenode = None # active node tracker
        self.activeedge = [] # active edge tracker
        self.activelength = 0 # active length tracker
        self.end = None #keep track of E
        self.nodedict = {} #dictionary to store the nodes
        self.leaves = 0 # leave counter
        self.inner_nodes = 0 # inner node counter
        self.suffixlink = [] # suffix link tracker
        self.splitpos = [] # split position 

    def increase_remainder(self):
        self.remainder = self.remainder + 1

    def decrease_remainder(self):
        self.remainder = self.remainder - 1

    def change_activenode(self,num):
        self.activenode = num
    
    def update_end(self, index):
        self.end = index
    
    def add_character_node(self, character, index):
        if character in self.nodedict:
            temp = self.nodedict[character]
            temp.append(index)
            self.nodedict[character] = temp
        else:
            self.nodedict[character] = [index]

            self.leaves = self.leaves + 1

    def show_all_nodes(self, T):
        for nodes in self.nodedict:
            index_list = self.nodedict[nodes]
            
            for node in index_list:
                print(T[node:self.end + 1])

    def change_activeedge(self, edge):
        self.activeedge.append(edge)

    def ret(self, leaveslist, innerleaveslist):
        leaveslist[6] = leaveslist[6] - 1
        innerleaveslist[6] = innerleaveslist[6] - 1
    
    def change_activelength(self):
        self.activelength = self.activelength + 1
    
    def decrease_activelength(self):
        self.activelength = self.activelength - 1

    def split(self,character, index):
        self.inner_nodes = self.inner_nodes + 1
        self.leaves = self.leaves + 1
        self.addsuffixlink(character)


    def clear_activeedge(self):
        self.activeedge = []

    def change_activenode(self, character):
        self.activenode = character

    def addsuffixlink(self, character):
        self.suffixlink.append(character)





def check_node_dict(node, character, index):
    if character in node.nodedict:
        if node.activenode is not None:
            node.activenode = None
            node.split(character, index)
            node.decrease_remainder()
        else:
            node.change_activeedge(character)            
            node.change_activelength()
            node.add_character_node(character, index)
            
        if character in node.suffixlink:
            node.change_activenode(character)
    else:
        if len(node.activeedge) != 0:
            for activeedge in node.activeedge:
                node.split(activeedge, index)
                node.decrease_remainder()
                node.decrease_activelength()
            
            node.clear_activeedge()

        node.add_character_node(character, index)
        node.decrease_remainder()





##################################################################

def build_suffixtree(T):
    """
    INPUT: Text T
    OUTPUT:
    - Your suffix tree ST
    - leaves: list of number of leaves after each step
    - inner_nodes: list of number of inner nodes after each step
    """
    ST = None # this should be your suffix tree
    leaves = [0]
    inner_nodes = [0]

    # checking the string conssits only ASCII characters > 37 except $
    for index in range(len(T)):
        if ord(T[index]) >= 37 or ord(T[len(T) - 1]) == 36:
            continue
        else:
            return 1

    d = node()
    leaveslist = []
    innerleaveslist = []

    for index, character in enumerate(T):
        d.increase_remainder()
        d.update_end(index)
        check_node_dict(d, character, index)
        d.show_all_nodes(T)
        leaveslist.append(d.leaves)
        innerleaveslist.append(d.inner_nodes)

    #for nodes in d.show_all_nodes(T):
    #    ST.append(nodes)

    d.ret(leaveslist, innerleaveslist)

    leaves = leaveslist
    inner_nodes = innerleaveslist
    ST = d

    return ST, leaves, inner_nodes

def get_text(args):
    if args.text is not None:
        return args.text
    with open(args.textfile, "r") as ftext:
        text = ftext.read()
    return text

def get_argument_parser():
    p = argparse.ArgumentParser(description="Ukonnen")
    txt = p.add_mutually_exclusive_group(required=True)
    txt.add_argument("-T", "--text",
        help="text to build suffix tree")
    txt.add_argument("-t", "--textfile",
        help="name of file containing text (will be read in one piece)")
    return p

def main(args):
    
    T = get_text(args)
    ST = test_suffixtree(T)

if __name__=="__main__":
    main(get_argument_parser().parse_args())
