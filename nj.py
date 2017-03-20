#!/usr/bin/python
# This is the python code for the Nussinov-Jacobson algorithm,
# the way we wrote it during class on April 4, 2016.

'''Goal: to find the optimal number of base pairs given an RNA sequence
Input: an RNA sequence
Ouput: the optimal number of base pair locations
'''

# define function isBasePair
def isBasePair(n1, n2) :
    if n1=='A' and n2=='U' :
        return 1
    elif n1=='U' and n2=='A' :
        return 1
    elif n1=='C' and n2=='G' :
        return 1
    elif n1=='G' and n2=='C' :
        return 1
    else :
        return 0

# define a function to determine the maximum k value
def getMaxK(matrix, i, j) :
    maxsum = 0
    # for k from i + 2 to j - 3
    for k in range(i+2, j-2) :
        # sum = the two matrix positions
        sum = matrix[i][k] + matrix[k+1][j]
        # if sum is more than maxsum then update maxsum
        if sum > maxsum :
            maxsum = sum 
    # return something
    return maxsum

# define a function to determine the optimal location of the base pairs
# this function is recursive
def findBasePairs(matrix, seq, r, c) :
    #base case: stop when row r is equal or larger than c
    if r >= c :
        return 
    # case C: when r, c spot is same as r+1, c spot, do the latter recursively
    if matrix[r][c] == matrix[r+1][c] :
        findBasePairs(matrix, seq, r+1, c)
    # case A: when r,c spot is same as to the left, do the left recursively
    elif matrix[r][c] == matrix[r][c-1] :
        findBasePairs(matrix, seq, r, c-1)
    # case B: when r,c spot is same as diagonally left plus the match
    elif matrix[r][c] == matrix[r+1][c-1] + isBasePair(seq[r],seq[c]) :
        # if the sequence at r basepairs with sequence at c
        if isBasePair(seq[r],seq[c]) :
            #print sequence at r (r) base pairs with sequence at c (c)
            print seq[r] + "(" + str(r) + ") base pairs with " + seq[c] + "(" + str(c) + ")"
        # recursively find the basepairs for cell r+1, c-1
        findBasePairs(matrix, seq, r+1, c-1)

    # else case D: 
    else :
        # for k from r+2 to c-2
        for k in range(r+2, c-1) :
            # let temp be matrix value at r,k plus at k+1,c
            temp = matrix[r][k] + matrix[k+1][c]
            # if temp is same as matrix at r,c
            if temp == matrix[r][c] :
                # then recursively find basepairs for these two parts.
                findBasePairs(matrix, seq, r, k)
                findBasePairs(matrix, seq, k+1, c)

# print the matrix
def printMatrix(matrix) :
    print "-----------------"
    print seq
    for row in matrix :
        print row
    print "-----------------"


# Step 1: initialize and read in the sequence
infile = open("sequence.txt", "r") 
line1 = infile.readline()
print "line1 = ", line1
seq = '' 
for line in infile :
    line = line.replace("\n", "") 
    seq = seq + line
print "sequence = ", seq
# let N be the length of the sequence
N = len(seq)
print "length = ", N
print "-----------------"

# Step 2: initialize an NxN matrix with zeros
matrix = [[0 for i in range(N)] for j in range(N)]

# Show the matrix
printMatrix(matrix)

# Step 3: fill in the matrix diagonal by diagonal
# for each diagonal from 1 to N-1
for diagonal in range(1, N) :
    # for each row from 0 to N - diagonal
    for row in range(N - diagonal) :
        # col is row + diagonal
        col = row + diagonal
        # w is what isBasePair gives for seq at row and seq at col
        w = isBasePair(seq[row], seq[col])
        # a =
        a = matrix[row][col-1]
        # b =
        b = matrix[row+1][col-1]
        # c =
        c = matrix[row+1][col]
        # d = getMaxK for this matrix at cell row, col
        d = getMaxK(matrix, row, col)
        # matrix at row, col equals the max of a, b, c, d
        matrix[row][col] = max(a, b+w, c, d)
# Show the matrix
printMatrix(matrix)

# use findBasePairs on the matrix, seq, starting at cell 0, N-1
# to print where the base pairs are
findBasePairs(matrix, seq, 0, N-1)

#-----------------------------------------
