import csv
import pickle
import numpy as np

print("Reading from GSE59793-GPL18998_series_matrix.txt")

f = open("GSE59793-GPL18998_series_matrix.txt","r")

entire = f.read()
rows = entire.split("\n")
data = []
for row in rows:
    data.append(row.split("\t"))

#for i in range(len(data)):
#    print("{0}: {1}".format(i,data[i][0]))

labels = data[27]
labels[0] = "Gene Name"

'''
Genes start at row 66. 66 - 198 are some mouse targeted genes (?)
199 - 5453 are some just genes. Just flat out genes.
Column entries are different samples at different time periods.

Let me try put it all in a numpy array.

So the array
Dim - What's there - |n|
0 - gene - ~5000
1 - sample - 2
2 - Channel - 2
3 - Stimulated - 2
4 - time - 10

'''


dtypes = []
for label in labels:
    if label == "Gene Name":
        dtypes += [(label, np.unicode_, 16)]
    else:
        dtypes += [(label, np.float)]
dtypes = np.dtype(dtypes)
#print(dtypes)
#raise Exception('')

X = np.empty(5388, dtype = dtypes)

l = len(labels)

for i in range(66, 5454):
    X[i-66][0] = data[i][0]
    for j in range(1, l):
        if data[i][j] == '':
            X[i-66][j] = 0
        else:
            X[i-66][j] = data[i][j]

with open('genedata.npy', 'wb') as f:
    np.save(f, X)
