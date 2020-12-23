
import numpy as np
import matplotlib.pyplot as plt

import csv

print("Performing discrimination")


#haha hoho lol
def haha(X, a, b):
    out = [X['"Replicate '+a+', Heavy channel SILAC, time 0"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 0.5h"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 1h"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 2h"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 3h"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 4h"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 5h"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 6h"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 9h"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 12h"'], X['"Replicate '+a+', Heavy channel SILAC, '+b+'-stimulus, time 24h"']]
    return out

def hoho(X, a, b):
    out = [X['"Replicate '+a+', Medium channel SILAC, time 0"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 0.5h"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 1h"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 2h"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 3h"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 4h"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 5h"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 6h"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 9h"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 12h"'], X['"Replicate '+a+', Medium channel SILAC, '+b+'-stimulus, time 24h"']]
    return out

def meanCurve(curveList):
    outList = np.empty((curveList.shape[1], curveList.shape[2]))
    for i in range(curveList.shape[1]):
        for j in range(curveList.shape[2]):
            entry = np.mean(curveList[:,i,j])
            outList[i,j] = entry
    return outList

def computeABC(curveA, curveB, xAxis):
    s = 0
    for i in range(len(xAxis)-1):
        ys = curveA[i] - curveB[i]
        yf = curveA[i+1] - curveB[i+1]
        h = (ys + yf)/2
        w = xAxis[i+1] - xAxis[i]
        s+= h*w
    s = np.abs(s)
    return s


def computeIndependences(meanCurves, timeSigs):
    independences = np.empty(meanCurves.shape[0])
    xRange = timeSigs[-1] - timeSigs[0]
    for i in range(meanCurves.shape[0]):
        minDist = 10 ** 8 #infinity would be great
        for j in range(meanCurves.shape[0]):
            if not i == j:
                d = computeABC(meanCurves[i], meanCurves[j], timeSigs) / xRange
                if d < minDist:
                    minDist = d
        independences[i] = minDist
    return independences

def computeVars(meanCurves, expressions):
    variances = np.empty(meanCurves.shape[0])
    for i in range(meanCurves.shape[0]):
        vS = 0
        particularMC = meanCurves[i]
        particularExpression = expressions[i]
        for timestep in range(particularMC.shape[0]):
            v = 0
            m = particularMC[timestep]
            e = particularExpression[:,timestep]
            for entry in e:
                v += (m - entry)**2
            v = np.sqrt(v / e.shape[0]) #variance at particular time-step
            vS += v
        vS = vS / particularMC.shape[0] #Mean of variances over all time-steps
        variances[i] = vS
    return variances

def discriminate(names, expressions, timeSigs):
    weightList = np.empty(expressions.shape[3])
#    print(expressions.shape) #Dims: LPS/MOCK (TYPE); Sample (R1 or R@); time; gene.
    meanCurves = np.empty((expressions.shape[0], expressions.shape[2], expressions.shape[3]))
    for i in range(expressions.shape[0]):
        MC = meanCurve(expressions[i])
        meanCurves[i] = MC
#    print(meanCurves.shape) #DIMS: mean of type, time, gene
#    print(weightList.shape)
    for i in range(expressions.shape[3]): #looping over all genes
        timeProfile = expressions[:,:,:,i] #sample type; sample; time
        meanProfile = meanCurves[:,:,i] #sample type; time
        independences = computeIndependences(meanProfile, timeSigs)
#        print("ID:{0}".format(i))
#        print(independences)
        variances = computeVars(meanProfile, timeProfile)
#        print(variances)
        metrics = np.empty(independences.shape[0])
        for j in range(metrics.shape[0]):
            metrics[j] = independences[j] / variances[j]
#        print(metrics)
        weight = np.mean(metrics)
#        print(weight)
#        print(meanProfile)
#        print("===")
        weightList[i] = weight
    return weightList

def plotshit(i):
    plt.title("Gene: {0}, rank: {1}, weight: {2}".format(geneNames[sortedArgs[i]], i, sortW[i]))
    LPSplot1 = [LPSR1[x][sortedArgs[i]] for x in range(11)]
    LPSplot2 = [LPSR2[x][sortedArgs[i]] for x in range(11)]
    MOCKplot1 = [MCKR1[x][sortedArgs[i]] for x in range(11)]
    MOCKplot2 = [MCKR2[x][sortedArgs[i]] for x in range(11)]

    plt.plot([0, 0.5, 1, 2, 3, 4, 5, 6, 9, 12, 24], LPSplot1,'r')
    plt.plot([0, 0.5, 1, 2, 3, 4, 5, 6, 9, 12, 24], LPSplot2,'r')
    plt.plot([0, 0.5, 1, 2, 3, 4, 5, 6, 9, 12, 24], MOCKplot1,'b')
    plt.plot([0, 0.5, 1, 2, 3, 4, 5, 6, 9, 12, 24], MOCKplot2,'b')
    
    plt.savefig('trajectory'+str(i).zfill(4)+'.png')
    plt.clf()


X = np.load(open('genedata.npy', 'rb'))
geneNames = X['Gene Name']

R1T0 = X['"Replicate 1, Heavy channel SILAC, time 0"']
R2T0 = X['"Replicate 2, Heavy channel SILAC, time 0"']

LPSR1 = haha(X, '1', 'LPS') #Heavy
LPSR2 = haha(X, '2', 'LPS') 
MCKR1 = haha(X, '1', 'Mock')
MCKR2 = haha(X, '2', 'Mock')

LPSR1H = hoho(X, '1', 'LPS') #Medium
LPSR2H = hoho(X, '2', 'LPS') 
MCKR1H = hoho(X, '1', 'Mock')
MCKR2H = hoho(X, '2', 'Mock')

timeSigs = [0, 0.5, 1, 2, 3, 4, 5, 6, 9, 12, 24]

weights = discriminate(geneNames,np.array([[LPSR1, LPSR2],[MCKR1, MCKR2]]), timeSigs)

#legitimateIndexes = np.argwhere(~np.isnan(weights)).T[0]
#legitimateNames = geneNames[legitimateIndexes]
#legitimateWeights = weights[legitimateIndexes]
#sortedArgs = np.argsort(legitimateWeights)[::-1]
#sortedWeights = legitimateWeights[sortedArgs]
#sortedNames = legitimateNames[sortedArgs]

#print(weights[:10])
args = np.array(range(weights.shape[0]))
#print(args[:10])
nonNan = np.argwhere(~np.isnan(weights)).T[0]
legitWeights = weights[nonNan]
legitArgs = args[nonNan]
#print(legitWeights[:10])
#print(legitArgs[:10])
sort = np.argsort(legitWeights)[::-1]
sortW = legitWeights[sort]
sortA = legitArgs[sort]
#print(sortW[:10])
#print(sortA[:10])
sortedArgs = sortA

LPSR1 = np.array(LPSR1)[:,sortedArgs]
LPSR2 = np.array(LPSR2)[:,sortedArgs]
MCKR1 = np.array(MCKR1)[:,sortedArgs]
MCKR2 = np.array(MCKR2)[:,sortedArgs]

LPSR1 = np.expand_dims(LPSR1, axis=0)
LPSR2 = np.expand_dims(LPSR2, axis=0)
MCKR1 = np.expand_dims(MCKR1, axis=0)
MCKR2 = np.expand_dims(MCKR2, axis=0)

LPSR1H = np.array(LPSR1H)[:,sortedArgs]
LPSR2H = np.array(LPSR2H)[:,sortedArgs]
MCKR1H = np.array(MCKR1H)[:,sortedArgs]
MCKR2H = np.array(MCKR2H)[:,sortedArgs]

LPSR1H = np.expand_dims(LPSR1H, axis=0)
LPSR2H = np.expand_dims(LPSR2H, axis=0)
MCKR1H = np.expand_dims(MCKR1H, axis=0)
MCKR2H = np.expand_dims(MCKR2H, axis=0)

sortedNames = geneNames[sortedArgs]
ready_data = np.empty((4, 2, LPSR1.shape[1], LPSR1.shape[2]),dtype = LPSR1.dtype) #Dims: 0-sample, type, 1-time, 2-gene

ready_data[0,0,:,:] = LPSR1
ready_data[1,0,:,:] = LPSR2
ready_data[2,0,:,:] = MCKR1
ready_data[3,0,:,:] = MCKR2
ready_data[0,1,:,:] = LPSR1H
ready_data[1,1,:,:] = LPSR2H
ready_data[2,1,:,:] = MCKR1H
ready_data[3,1,:,:] = MCKR2H


np.save(open('sorted_genes.npy', 'wb'), ready_data)
np.save(open('sorted_names.npy', 'wb'), sortedNames)
sortedX = X[sortedArgs]
