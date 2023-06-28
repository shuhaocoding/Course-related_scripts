# Algorithms for DNA Sequencing W1
import os
os.chdir('/Users/shuhao/Desktop')
def naive(p, t):
    a=0;b=len(t) - len(p) + 1
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):
            a+=1 # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences,a,b
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def naiveRC(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
        if reverseComplement(p) != p:
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i + j] != reverseComplement(p)[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
    return occurrences

genome = readGenome('lambda_virus.fa')
min(naiveRC('AGTCGA',genome))

def naivemis2(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        count = 0
        for j in range(len(p)):
            if t[i+j] != p[j]:  # compare characters
                count += 1
                if count > 2:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

min(naivemis2('AGGAGGTT', genome))

seq, qual = readFastq('ERR037900_1.first1000.fastq')
scores = []

for i in range(len(qual[0])):
    s = 0
    for sq in qual:
        s += ord(sq[i])
    scores.append(s)

print(scores.index(min(scores)))

# Algorithms for DNA Sequencing W2
def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    a=0;b=0
    occurrences = []
    while i < len(t) - len(p) + 1:
        b+=1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            a+=1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences,a,b
def naive(p, t):
    a=0;b=len(t) - len(p) + 1
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):
            a+=1 # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences,a,b
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
t = readGenome('chr1.GRCh38.excerpt.fasta')

occurrences,a,b = naive('GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG', t)

p = 'GGCGCGGTGGCTCACGCCTGTAAT'
p_bm = BoyerMoore(p)
occurrences,a,b = boyer_moore(p, p_bm, t)

p = 'GGCGCGGTGGCTCACGCCTGTAAT'

c1 = int(len(p)//3)
c2 = int(len(p)//3 *2)
p1=p[:c1]
p2=p[c1:c2]
p3=p[c2:]

t_index = Index(t,8)
i1 = t_index.query(p1)
i2 = [i-8 for i in t_index.query(p2)]
i3 = [i-16 for i in t_index.query(p3)]

s = set(i1)|set(i2)|set(i3)
print(len(s))

print(len(i1)+len(i2)+len(i3))

p = 'GGCGCGGTGGCTCACGCCTGTAAT'
t_subindex = SubseqIndex(t,8,3)
i1 = t_subindex.query(p[0:])
i2 = [i-1 for i in t_subindex.query(p[1:])]
i3 = [i-2 for i in t_subindex.query(p[2:])]
print(len(i1)+len(i2)+len(i3))

# Algorithms for DNA Sequencing W3
import os
os.chdir('/Users/shuhao/Desktop')
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]
def read_FAST_A(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
def editDistanceDP(A, B):
    m, n = len(A), len(B)
    dp = [[0 for j in range(n+1)] for i in range(m+1)]
    for i in range(m+1): dp[i][0] = i
    for j in range(n+1): dp[0][j] = j
    for i in range(1, m+1):
        for j in range(1, n+1):
            dp[i][j] = min(
                dp[i-1][j-1] + int(A[i-1] != B[j-1]),
                dp[i-1][j] + 1,
                dp[i][j-1] + 1,
            )
    return dp[m][n]
def editDistanceApproximate(P, T):
    m, n = len(P), len(T)
    dp = [[0 for j in range(n + 1)] for i in range(m + 1)]
    for i in range(m + 1): dp[i][0] = i  # init first column by distance from empty string

    #   the first row is all 0s unlike edit distance, since there is no bias toward alignment
    #   of P in T from the beginning of both P and T, (ie. P can start at any index in T)

    #
    #   for j in range(n+1): dp[0][j] = j # DELETED!!!
    #

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            dp[i][j] = min(
                dp[i - 1][j - 1] + int(P[i - 1] != T[j - 1]),
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
            )
    return min(dp[m])
P ='GCGTATGC'
T ='TATTGGCTATACGGTT'
print(editDistanceApproximate(P, T))
T = readGenome('chr1.GRCh38.excerpt.fasta')
P = 'GATTTACCAGATTGAG'
print(editDistanceApproximate(P, T))

#W4
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

import itertools

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest
i = ['CCT','CTT','TGC','TGG', 'GAT','ATT']
print(scs(i))