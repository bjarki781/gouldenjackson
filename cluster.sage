from sage.combinat.species.generating_series import *

# note that here the variable T is always the set of subwords/subsequences
# that we want to find the distribution of.
# m is the size of T's alphabet 
# w is an index into T which signifies the word we are focussing on,
# the winning word in penney's game

# test to see if a set of words are reduced
def is_reduced(T):
    for t in T:
        for s in T - {t}:
            if t.find(s) != -1:
                return false
    return true

# the weighted adjacency matrix of the words, how they can cluster together
def wam(T):
    if not is_reduced(set(T)):
        print("Patterns given aren't reduced")
        return []

    list = []
    us = []

    line = [0]
    for i in range(1, len(T)+1):
        us.append(var('u' + str(i)))
        line.append(us[i-1] * x^len(T[i-1]))
    list.append(line)

    for i in range(1, len(T)+1):
        line = [0]
        for j in range(1, len(T)+1):
            line.append(us[j-1]*psi_x(T[i-1], T[j-1]))
        list.append(line)
    return matrix(list)

# the matching operation used when calculation correlation between
# two strings
def match(p, q, i):
    n = len(p)-i
    if p[i:] == q[:n]:
        return true
    return false

# defined as Omega in my text and in Guibas and Odzlyko PQ
def correlation(p, q):
    s = ""
    for i in range(1, len(p)):
        if match(p, q, i):
            s += "1"
        else:
            s += "0"
    return s

# my altered omega function, so we can cluster words of differing lengths
def psi_x(p, q):
    r = 0
    om = correlation(p,q)
    if len(p) < len(q):
        zeros = zero_string(len(p)-len(q))
        om = zeros + om
    elif len(p) > len(q):
        om = om[len(p) - len(q):]

    for i in range(0, len(om)):
        if om[i] == "1":
            r += x^(i+1)
    return r

# utility function for psi_x
def zero_string(n):
    s = ""
    for i in range(0, n):
        s += "0"
    return s

# the normal cluster generating function
def clustergf(T):
    B = wam(T)
    n = B.nrows()
    R = ~(Matrix.identity(n)-B)

    return sum(R.row(0)[1:])

# our altered cluster generating function, C_w in my text
def penneysclustergf(T, s):
    B = wam(T)
    n = B.nrows()

    ui = var('u' + str(s+1))

    B = insert_row(B, n, zero_vector(n))
    B = insert_column(B, n, B.column(s+1)/ui)

    # n+1 since we added a column
    R = ~(Matrix.identity(n+1)-B)

    return sum(R.row(0)[1:])

# utility operations 
def insert_row(M,k,row):
    return matrix(M.rows()[:k]+[row]+M.rows()[k:])
def insert_column(M,k,column):
    return transpose(matrix(M.columns()[:k]+[column]+M.columns()[k:]))

# the normal distribution of subwords
# m is the size of the alphabet
def dist(T, m):
    Cw = penneysclustergf(T, s).substitute(submap)
    C = clustergf(T).substitute(submap)

    return 1/(1 - (m*x + C)).full_simplify()

# distribution of subwords that end in an unmarked word
# m is the size of the alphabet
def penneysdist(T, m, w):
    submap = {}
    for i in range(0, len(T)):
        submap.update({var('u' + str(i+1)):-1})
    Cw = penneysclustergf(T, w).substitute(submap)
    C = clustergf(T).substitute(submap)

    return (Cw/(1 -(m*x + C))).full_simplify()


# calculate the odds of T[w] winning penney's game
# with an m-sided dice,
# n is the amount of the terms expanded, higher n => more accuracy
def sequence_odds(T, m, w, n):
    ce = penneysdist(T, m, w).series(x, n).list()
    pc = []
    for i in range(0, len(ce)):
        pc.append(ce[i]/m^i)

    return pc

def all_odds(T, m=2):
    r = []
    # a guesstimate of terms needed to get an acceptably accurate answer
    n = max(map(len, T))*5^m
    for w in range(0, len(T)):
        r.append(1 + sum(sequence_odds(T, m, w, n)).n())
    return r

