import numpy as np

Base2Index_DNA = {'A':0, 'C':1, 'T':2, 'G':3}
Base2Index_Bin = {'0':0, '1':1}
def cgk_embedding(s, N, R, seq_type=0, padding=True):

    #pdb.set_trace()

    if seq_type==1:
        s = s.upper()

    #L_R = 3*N*4
    #R = np.random.randint(0,2,L_R)
    
    i = 0
    s1 = ""
    for j in range(3*N):
        if i<len(s):
            c = s[i]
            s1 = s1 + c
            if seq_type==0:
                i = i + R[(j-1)*2+Base2Index_Bin[c]]
            else:
                i = i + R[(j-1)*4+Base2Index_DNA[c]]
        elif padding==True:
            s1 = s1 + "N"
        else:
            break
    return s1

def hamming(a,b):
    return sum([1 for i in range(len(a)) if a[i]!=b[i]])

'''
proj a --> a1, and b --> b1
to use ham(a1, b1) to bound edit(a,b)
'''
def proj_hamming_dist_bkp(a, b, seq_type, normalize=True):
    #pdb.set_trace()
    N_iter = 10
    N = max(len(a), len(b))
    min_dist = 3*N+1
    for i in range(N_iter):

        if seq_type==0:
            L_R = 3*N*2
        else:
            L_R = 3*N*4
        R = np.random.randint(0,2,L_R)

        a1 = cgk_embedding(a, N, R, seq_type)
        b1 = cgk_embedding(b, N, R, seq_type)
        h_a1_b1 = hamming(a1, b1)#a1,b1 of 3N lengths
        #print(h_a1_b1)
        if h_a1_b1<min_dist:
            min_dist = h_a1_b1

    if normalize==True:
        min_dist = float(min_dist)/(3*N)
    else:
        min_dist = float(min_dist)
    #print(min_dist)
    return min_dist

'''
proj a --> a1, and b --> b1
to use ham(a1, b1) to bound edit(a,b)

R_list: can be None or not
'''
def proj_hamming_dist(a, b, seq_type, normalize=True,
                      R_list=None):
    N = max(len(a), len(b))
    #pdb.set_trace()
    if R_list is None:
        N_iter = 10

        if seq_type==0:
            L_R = 3*N*2
        else:
            L_R = 3*N*4

        R_list = []
        for i in range(N_iter):
            R = np.random.randint(0,2,L_R)
            R_list.append(R)

    min_dist = 3*N+1
    for R in R_list:
    
        a1 = cgk_embedding(a, N, R, seq_type)
        b1 = cgk_embedding(b, N, R, seq_type)
        h_a1_b1 = hamming(a1, b1)#a1,b1 of 3N lengths
        #print(h_a1_b1)
        if h_a1_b1<min_dist:
            min_dist = h_a1_b1

    if normalize==True:
        min_dist = float(min_dist)/(3*N)
    else:
        min_dist = float(min_dist)
    #print(min_dist)
    return min_dist

np.random.seed(0)

a = "01"
b = "0101"

N = max(len(a),len(b))
L_R = 3*N*2

R_list = []
N_iter = 1
for i in range(N_iter):
    R = np.random.randint(0,2,L_R)
    R_list.append(R)

print(proj_hamming_dist(a,b,0,False, R_list))
