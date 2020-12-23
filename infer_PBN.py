import numpy as np
import pickle
import gym
import networkx as nx
import matplotlib.pyplot as plt
import math
print("")
print("Inferring PBN")

def integerize(X):
    n_in = X.shape[0]
    s = [0,0]
    for i in range(n_in):
        if X[n_in-i-1] == 0.5:
            s[0] += 0 * 2**i
            s[1] += 1 * 2**i
        else:
            s[0] += int(X[n_in-i-1]) * 2**i
            s[1] += int(X[n_in-i-1]) * 2**i
    return s

def booleanize(X, l):
    output = np.zeros((1,l), dtype=bool)
    for i in range(l):
        h = 2**(l-i-1)
        if X >= h:
            X -= h
            output[0,i] = 1
    return output


def lin_reg(X,Y):
    n_exp, n_in = X.shape
    ones = np.ones((n_exp, 1), dtype=bool)
    X = np.append(X, ones, axis=1).astype(int)
    R = np.matmul(X.T, X)
    Rp = np.linalg.pinv(R)
    C = np.matmul(X.T,Y)
    A = np.matmul(Rp, C)
    return A

def construct_inputs(length):
    output = np.empty((2**length, length),dtype=bool)
    for i in range(2**length):
        b = booleanize(i, length)
        output[i,:] = b
    return output

def construct_F_A(A):
    n_in = A.shape[0]-1
    X = construct_inputs(n_in)
    n_samples, n_in = X.shape
    ones = np.ones((n_samples, 1), dtype=bool)
    X = np.append(X, ones, axis=1)

    output = np.matmul(X,A)
    return output

#Finding opposite corner of a hypercube
def find_opposite(index, dim):
    bool_expr = ~booleanize(index, dim)[0,:]
    opp = integerize(bool_expr)[0]
    return opp

def normalize(F):
    dim = int(math.log2(F.shape[0]))
    neg_flag = F<0
    pos_flag = F>1
    if (neg_flag == True).any():
        indexes = np.arange(F.shape[0])
        neg_indexes = indexes[neg_flag]
        for index in neg_indexes:
            opp = find_opposite(index, dim)
            err = -F[index]
            F[index] += err
            F[opp] -= err
    if (pos_flag == True).any():
        indexes = np.arange(F.shape[0])
        pos_indexes = indexes[pos_flag]
        for index in pos_indexes:
            opp = find_opposite(index, dim)
            err = 1-F[index]
            F[index] += err
            F[opp] -= err
    return F

def construct_F(X,Y):
    n_exp, n_in = X.shape
    F = np.zeros((2**n_in), dtype=int)
    sums = np.zeros((2**n_in), dtype=int)
    for exp_i in range(n_exp):
        x = X[exp_i, :]
        index = integerize(x)
        sums[index] += 1
        if Y[exp_i] == 1:
            F[index] += 1

    if (sums==0).any():
        F = np.divide(F, sums)
        A = lin_reg(X,Y)
        F_lin = construct_F_A(A)
        F_norm = normalize(F_lin[:,0])
#        print(F)
#        print(F_lin[:,0])
#        print(F_norm)
        F = F_norm
#        print("---")
#        raise RuntimeWarning("This needs fixing since linalg may fuck up")
    else:
        F = np.divide(F, sums)
    return F

def find_top_input_gene(experience, mask): #return: value, F, mask
    n_genes = mask.shape[0]
    n_exp = experience.shape[0]
    n_in = np.sum(mask)
    X = np.empty((n_exp, n_genes), dtype=bool)
    Y = np.empty((n_exp,1), dtype=bool)
    for i in range(n_exp):
        Y[i,0] = experience[i][1]
        x = experience[i][0]
        X[i] = x

    potential_inputs = np.arange(n_genes)[np.invert(mask)]

    value_max = 0
    F_max = None
    mask_max = None
    delta_max = None

    for pot_i in potential_inputs:

        consider_mask = np.copy(mask)
        consider_mask[pot_i] = True

        X_consider = X[:,consider_mask]
        F_new = construct_F(X_consider, Y)
        value = np.max(F_new) - np.min(F_new)
#        if value > 1:
#            print(F_new)
#            f_norm = normalize(F_new)
#            print(f_norm)
#            raise Exception('')
        if value > value_max:
            value_max = value
            F_max = F_new
            mask_max = consider_mask
#        print("Function: {0}".format(np.around(F_new,3)))
#        print("Convergence: {0}".format(np.max(F_new)))
#        print("Divergence: {0}".format(np.min(F_new)))
#        print("Value: {0}".format(value))
#        print("------")
#    print("Max: {0}".format(value_max))
#    print("=====")
    return value_max, F_max, mask_max


genedata = np.load(open("quantized_sorted_genes.npy", "rb"))
genenames = np.load(open("sorted_names.npy", "rb"))

n_genes = 20


genedata = genedata[:,:,:n_genes]
genenames = genenames[:n_genes]

n_samples, n_time, n_genes = genedata.shape

#Tuplerizing

n_exp = n_samples * (n_time-1)

experiences = np.empty((n_exp, n_genes), dtype=object)

for i in range(n_genes):

    for t in range(n_time-1):
        for s in range(n_samples):
            target_expr = genedata[s,t,i]
            input_expr = genedata[s,t+1, :]
            tup = (input_expr, target_expr)
            experiences[t*(n_samples) + s,i] = tup

#Woohoo. Now PBN inference.


PBN_data = np.empty((n_genes), dtype=object)

value_min = 10**-8

PBN_data = []

for gene_i in range(n_genes):
    input_mask = np.zeros((n_genes), dtype=bool)
    e = experiences[:,gene_i]
    
    value_old = 0

    value_new, F, input_mask = find_top_input_gene(e, input_mask)
    while value_new <=1 and (value_new - value_old) > value_min:
        print("Gene {1}: Adding  gene with value {0}".format(value_new, gene_i))
        value_old = value_new
#        print(value_old)
        input_mask_old = input_mask
        F_old = F
        value_new, F, input_mask= find_top_input_gene(e, input_mask)
#    print("Done")
    PBN_data += [(input_mask_old, np.around(F_old,4))]
#    print(value_old)
#    print(input_mask_old.astype(int))
#    print(np.around(F_old,3))

#print(PBN_data)
#print(len(PBN_data))
environment = gym.make('gym_PBN:PBN-v0', PBN_data = PBN_data)

pickle.dump(environment, open("PBN_jovanovic.pkl", "wb"))
