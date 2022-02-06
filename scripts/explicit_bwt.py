#!/usr/bin/env python3 



def rotations(T):
    """ This method returns a list of rotations for a string T"""
    ## Take twice the string then slide a window of len(string) 
    TT = T*2
    return [ TT[i:i+len(T)] for i in range(0, len(T))]

def get_bwm(T):
    """ Sort rotations lexicographically"""

    return sorted(rotations(T))


def BWT(T):
    """Gives the BW transformed string BWT(T)"""
    ## Map function onto all elements of bwm list
    ## basically only takes the last element 

    return ''.join(map(lambda x : x[-1], get_bwm(T)))


def print_rotations(matrix):
    """Prints the rotations and BW matrix"""
   
    for rotation in matrix:
        for i in range(len(rotation)):
            print(rotation[i]," ", end='')
        print("")

    print('\n----------------------')

def print_BWT(bwt):
    """Prints transformed T string"""
    for char in bwt:
        print(char,' ', end= ' ')
    print('\n----------------------')


T = "$acaacg"

print("\nThe Rotations:")
print_rotations(rotations(T))


print("The BWM:")
print_rotations(get_bwm(T))


transformed = BWT(T)
print("The transformed String:")
print_BWT(BWT(T))

