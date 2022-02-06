#!/usr/bin/env python3 


import sys
import os



args = sys.argv
if len(args) == 1:
    print('Please Specify filename containing sequence\nUsage: ./script.py <path/to/file>')
    sys.exit()

try:
    f = open(args[1], 'rt')
except:
    print('File Not Found!')
    sys.exit()


for line in f:
    T = line.strip()
    break


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



def print_BWT(bwt):
    """Prints transformed T string"""
    for char in bwt:
        print(char,' ', end= ' ')
    print('\n-------------------------')


transformed = BWT(T)
print("\nThe transformed String:")
print('-------------------------')
print_BWT(transformed)

