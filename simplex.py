#!/usr/bin/python3
# Author: Suzanna Sia

# Standard imports
#import random
import numpy as np
import pdb
import math
import os, sys

# argparser
#import argparse
#from distutils.util import str2bool
#argparser = argparser.ArgumentParser()
#argparser.add_argument('--x', type=float, default=0)

# Custom imports

def construct_tableau(A, b, c):

    b = np.append(0, b)
    b = b.reshape(1, b.shape[0])
    zcol = np.zeros((A.shape[0]+1))
    zcol[0] = 1
    zcol = zcol.reshape(1, zcol.shape[0])
    
    pret = np.vstack((-c, A))
    prett = np.hstack((zcol.T, pret, b.T))

    return prett

def row_echelon(prett, pivots):
    #for (i, col) in enumerate(basis_columns):
    for (row, col) in pivots:
        #row = i+1 # index + 1 because we don't touch the first row
        prett[row] = prett[row]/prett[row, col]
        
        for irow in range(prett.shape[0]):
            if irow!=row:
                mul = prett[irow, col]/prett[row,col]
                prett[irow] -= prett[row]*mul

    return prett

def find_pivot(prett):

    col = np.argmax(prett[0][1:-1])+1
    ratios = prett[1:, -1]/prett[1:, col]
    ignorec = np.where(ratios <=0)[0]
    ratios[ignorec] = np.inf

    row = np.argmin(ratios) + 1 # for first row
    return (row, col)

def get_bfs(prett):

    columns = np.where(prett[0]==0)[0]
    bfs = np.zeros(prett.shape[1]-1)
    for col in columns:
        row = np.argmax(prett[:,col])
        val = prett[row, -1]
        bfs[col-1] = val
    
    return bfs

if __name__ == "__main__":
# 1. construct pretableau
# 2. row_echelon on basis columns
# 3. identify new basis column and pivot row position
    prob=3

    A = np.loadtxt(f'inputs/A{prob}.txt')
    b = np.loadtxt(f'inputs/b{prob}.txt')
    c = np.loadtxt(f'inputs/c{prob}.txt')
    pivots = np.loadtxt(f'inputs/pivots{prob}.txt')
    pivots = [(i+1, int(p)) for i, p in enumerate(pivots)]
   
    prett = construct_tableau(A, b, c)
    prett = row_echelon(prett, pivots)

    while len(np.where(prett[0][1:-1]<=0)[0])!=A.shape[1]:
        new_piv = find_pivot(prett)
        prett = row_echelon(prett, [new_piv])

    bfs = get_bfs(prett)
    objv = prett[0, -1]

    print("bfs:", bfs)
    print("objv:", objv)

