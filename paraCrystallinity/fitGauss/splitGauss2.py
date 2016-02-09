import os
import numpy as np
import pylab as P
import random as R
import csv
import mpl_toolkits.mplot3d.axes3d as p3
import time as T
import copy
import sys
sys.path.insert(0, '/panfs/panasas1.hpc.dur.ac.uk/ghsk28/CustomPythonModules')
#import scipy.stats as stats
from lmfit.models import GaussianModel

def getFilesList(direc):
    fileList = os.listdir(direc)
#    files = []
#    for element in fileList:
#        files.append(element[7:])
    csvFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.csv'):
            csvFiles.append(str(fileName))
    return csvFiles

def readCSV(direc, csvFile):
    dhklData = []
    csvHandle = open(direc+csvFile)
    csvRead = csv.reader(csvHandle)
    for row in csvRead:
        dhklData.append(float(row[0]))
    csvHandle.close()
    return dhklData

def analyseComponents(gausses):
    components = gausses.components
    componentData = []
    for i in range(len(components)):
        componentData.append([])
        characterString = ''
        recordValue = False
        for character in str(components[i]):
            if (recordValue == True) and (character != ' ') and (character != ']') and (character != ','):
                characterString += character
            if character == '[':
                recordValue = True
            elif character == ',':
                recordValue = True
                componentData[-1].append(float(characterString))
                characterString = ''
            elif character == ']':
                componentData[-1].append(float(characterString))
        componentData[-1].append(gausses.pi[i])
    return componentData


def gaussianForceMag(x, mu, sigma, desiredMagnitude):
    gaussian = (1/(np.sqrt(2*np.pi*sigma**2)))*np.exp(- (x - mu)**2/(2*sigma**2))
    scaleFactor = desiredMagnitude/np.max(gaussian)
    return gaussian*scaleFactor

def gaussian(x, mu, sigma, scaleFactor):
    gaussian = scaleFactor*(1/(np.sqrt(2*np.pi*sigma**2)))*np.exp(- (x - mu)**2/(2*sigma**2))
    return gaussian

def plotDist(dhklData, componentData, forceMagnitude, csvFileName):
    fileName = csvFileName[:-4]
    P.figure()
    n, bins, patches = P.hist(dhklData, bins=1000, normed=1)
    P.xlabel('Dhkl')
    P.ylabel('Occurence')
    P.xlim([0,30])
    P.savefig('./outputFiles/'+fileName+'.png')
#    dhklData.sort()

    for i in bins:
        if bins[i] > 6:
            break
    firstGaussMax = np.max(n[:i])
    secondGaussMax = np.max(n[i:])

    # print "---==EXPT===---"
    # print firstGaussMax
    # print secondGaussMax
    # print "---======---"

    xvals = np.arange(0,30,0.01)

    if forceMagnitude == True:
        firstGauss = gaussianForceMag(xvals, componentData[0][0], componentData[0][1], firstGaussMax)
        secondGauss = gaussianForceMag(xvals, componentData[1][0], componentData[1][1], secondGaussMax)
    else:
        firstGauss = gaussian(xvals, componentData[0][0], componentData[0][1], componentData[0][2])
        secondGauss = gaussian(xvals, componentData[1][0], componentData[1][1], componentData[1][2])

    # print "---==THEORY===---"
    # print np.max(firstGauss)
    # print np.max(secondGauss)
    # print "---======---"



#    print firstGauss

    print "Number of dhkls =", len(dhklData)
    print "Trapz of dhkl dist =", np.trapz(n)
    print firstGauss
    print "Trapz of firstGauss =", np.trapz(firstGauss, xvals)
 


#    firstGauss = stats.norm.pdf(bins, componentData[0][0], componentData[0][1])*componentData[0][2]
#    secondGauss = stats.norm.pdf(bins, componentData[1][0], componentData[1][1])*componentData[1][2]
    P.plot(xvals, firstGauss, linewidth=2.0)
    P.plot(xvals, secondGauss, linewidth=2.0)
    if (forceMagnitude != True):
        P.plot(np.arange(0,30,0.01), firstGauss+secondGauss, linewidth=2.0, c='c', linestyle='--')
#    P.plot(gaussian(componentData[1][0], componentData[1][1]))
    P.savefig('./outputFiles/'+fileName+'_gauss.png')

    if (forceMagnitude == True):
        print "WARNING! THE MAGNITUDE FORCER IS ACTIVE, WHICH WILL BREAK THE FOLLOWING SUBROUTINES! TURN IT OFF!"


    for i in range(len(xvals)):
        if (xvals[i] > componentData[0][0]): # If we're past the mean of the first gaussian
            if (secondGauss[i] > firstGauss[i]): # This is the crossing point to the amorphous region
                break

    intersectPoint = xvals[i]


    thiosInCrystal = []
    for dhklVal in dhklData:
        if dhklVal <= intersectPoint:
            thiosInCrystal.append(dhklVal)



    print "\n\n----==== DISTRIBUTION RESULTS ====----"
    print "There are", len(thiosInCrystal), "thiophenes with neighbours within the intersect point ("+str(intersectPoint), "angstroms), out of", len(dhklData), "thiophenes considered, representing a portion of", str(float(len(thiosInCrystal))/float(len(dhklData)))+"."


    crystalArray = np.array(thiosInCrystal)
    crystalMean = np.mean(crystalArray)
    crystalStandardDev = np.std(crystalArray)

    print "\n"

    print "The two Gaussian distributions are mixed in a ratio of crystalline:amorphous of", str(componentData[0][2])+":"+str(componentData[1][2])+"."


    print "\n"
    print "Using the crossing point of the two gaussians leads to the following data:"
    print "Mean dhkl Value =", crystalMean
    print "Standard Deviation =", crystalStandardDev
    print "Paracrystallinity (using s^{2} = d^{2} * g^{2}) =", crystalStandardDev/float(crystalMean)

    print "\n"
    print "Using the crystalline gaussian distribution only (more accurate!) we get:"
    print "Mean dhkl Value =", componentData[0][0]
    print "Standard Deviation =", componentData[0][1]
    g = componentData[0][1]/float(componentData[0][0])
    print "Paracrystallinity (using s^{2} = d^{2} * g^{2}) =", g

    return g

def writeCSV(xData, yData, name):
    '''Writes the data values in a CSV file'''
    filename = './outputFiles/'+name+'.csv'
    document = csv.writer(open(filename, 'w+'), delimiter = ',')
    for i in range(len(xData)):
        document.writerow([xData[i], yData[i]])
    print 'Data written to', filename


def plotPara(timestepNos, paraData, filename):
    P.figure()
    P.plot(timestepNos, paraData)
    P.xlabel('Timestep')
    P.ylabel('Paracrystallinity, g')
    P.savefig('./outputFiles/para_'+str(filename)+'.png')


if __name__ == '__main__':

    forceMagnitude = False # Renomarlise Gaussians for easier comparison of shape



    direc = './inputFiles/'
    csvFiles = getFilesList(direc)
    csvFiles.sort()
    print "Files to treat:", csvFiles

    print len(csvFiles), "csv files found..."


    paraFileName = ''

    timestepNos = []
    paraCrystallinity = []

    for csvFile in csvFiles:
        print "Reading in", str(csvFile)+'...'

        stepNumber = int(csvFile[-7:-4])
        currentFileName = csvFile[:-12]

        # Make the paracrystallinity plot if we're looking at a different file

        paraFileName = currentFileName

        dhklData = readCSV(direc, csvFile)

        # Use PyMix to figure out the gaussians

     

        dhklUnder6 = []
        dhklOver6 = []
        for dhklVal in dhklData:
            if dhklVal < 6:
                dhklUnder6.append(dhklVal)
            else:
                dhklOver6.append(dhklVal)

        print "Guess for first gaussian", np.mean(dhklUnder6), np.std(dhklUnder6)
        print "Guess for second gaussian", np.mean(dhklOver6), np.std(dhklOver6)


        xvals = np.arange(0,30,0.01)
        gauss1 = GaussianModel()
        pars = gauss1.guess(dhklData, x=xvals)
        out = gauss1.fit(dhklData, pars, x=xvals)

#        result = gauss1.fit(dhklData, x=xvals, mu=np.mean(dhklData), sigma=np.std(dhklData), scaleFactor=1)

        print out.fit_report(min_correl=0.25)

        raise SystemError('STOP')





        # Insert guesses for Normal Distributions and Mixture
        n1 = mixture.NormalDistribution(np.mean(dhklUnder6), np.std(dhklUnder6)) # (mu, sigma)
        n2 = mixture.NormalDistribution(np.mean(dhklOver6), np.std(dhklOver6)) # (mu, sigma)


        ##### Dynamically fit 2 gaussians
        data = mixture.DataSet()
        data.fromList(dhklData)

        m = mixture.MixtureModel(2, [0.6,0.4], [n1, n2]) # (number of Gaussians, relative weights, Gaussian objects)

        # Perform expectation maximisation
        m.EM(data, 100, 0.00001) # (data, maximum iterations, convergence tolerance)

        print m
        gaussModel = analyseComponents(m) # [[crystalMean, crystalSD, crystalProp], [amorphMean, amorphSD, amorphProp]]

        paraG = plotDist(dhklData, gaussModel, forceMagnitude, csvFile)

        if (paraG >= 0.3):
            #### Forcibly fit 2 gaussians
            print "\n\n--------------============================================---------------"
            print "ONLY ONE GAUSSIAN FOUND - MORPHOLOGY IS AMORPHOUS AT THIS TEMPERATURE"
            print "FORCING TWO GAUSSIANS BY SPLITTING DATA AROUND d = 6 ANGSTROMS"
            print "--------------============================================---------------\n\n"
            data1 = mixture.DataSet()
            data2 = mixture.DataSet()
            data1.fromList(dhklUnder6)
            data2.fromList(dhklOver6)
            m1 = mixture.MixtureModel(1, [1.0], [n1])
            m2 = mixture.MixtureModel(1, [1.0], [n2])

            m1.EM(data1, 100, 1E-5)
            m2.EM(data2, 100, 1E-5)

            print m1
            print m2

            gaussModel1 = analyseComponents(m1) # [[crystalMean, crystalSD, crystalProp], [amorphMean, amorphSD, amorphProp]]        
            gaussModel2 = analyseComponents(m2) # [[crystalMean, crystalSD, crystalProp], [amorphMean, amorphSD, amorphProp]]
            gaussModel = [gaussModel1[0], gaussModel2[0]]


        paraG = plotDist(dhklData, gaussModel, forceMagnitude, csvFile)

        timestepNos.append(stepNumber)
        paraCrystallinity.append(paraG)


        if ((len(timestepNos) != 0) and (currentFileName != paraFileName)) or (csvFile == csvFiles[-1]):
            print "Plotting graph"
            print "Timesteps =", timestepNos
            print "g =", paraCrystallinity
            paraFileName = currentFileName
            plotPara(timestepNos, paraCrystallinity, paraFileName)
            writeCSV(timestepNos, paraCrystallinity, paraFileName)
            timestepNos = []
            paraCrystallinity = []
            



