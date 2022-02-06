#!/usr/bin/env python3 



def rotations(T):
    """ This method returns a list of rotations for a string T """
    ## Take twice the string then slide a window of len(string) 
    TT = T*2
    return [ TT[i:i+len(T)] for i in range(0, len(T))]

def get_bwm(T):
    """ Sort rotations lexicographically"""

    return sorted(rotations(T))


def BWT(T):
    """Gives the BW transformed string BWT(T)"""
    
    return ''.join(map(lambda x: x[-1], get_bwm(T)))


def print_rotations(matrix):
    """Prints the rotations and BW matrix"""
   
    for rotation in matrix:
        for i in range(len(rotation)):
            print(rotation[i]," ", end='')
        print("")



my_list = rotations("$agac")

print_rotations(my_list)
