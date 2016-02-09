import os
import numpy as np
import pylab as P
import time as T
import csv
import subprocess
import copy


def getInpList(direc):
    fileList = os.listdir(direc)
    inpFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.inp'):
            inpFiles.append(str(fileName))
    return inpFiles


def determineMonomersInSegment(path):
    inputHandle = open(path, 'r')
    inputLines = inputHandle.readlines()
    inputHandle.close()
    numberOfMonomers = 0
    recordAtoms = False
    for line in inputLines:
        if "xyz" in line:
            recordAtoms = True
            continue
        if (recordAtoms == True) and ('*' in line):
            recordAtoms = False
            break
        if recordAtoms == True:
            test = line.split(' ')
            if test[1] == 'S':
                numberOfMonomers += 1
    return numberOfMonomers


def obtainSegmentNumbers(name):
    segmentNumbers = []
    numberOfSegments = name.count("seg")
    test = name.split("seg")
    if (numberOfSegments == 1):
        for element in test:
            if (len(element) != 0):
                segmentNumbers.append(int(element))
    elif (numberOfSegments == 2):
        segmentNumbers = []
        for element in test:
            if (len(element) != 0):
                if (element[-1] == "_"):
                    segmentNumbers.append(int(element[:-1]))
                else:
                    segmentNumbers.append(int(element))
    return segmentNumbers


if __name__ == "__main__":

    validMorphs = os.listdir('./orcaInputs')
    validMorphs = sorted(validMorphs)

    completedInputFiles = 0
    totalInputFiles = 0
    t0 = T.time()

    for morphology in validMorphs:
        fullPath = './orcaInputs/'+morphology+'/'
        validInpsList = getInpList(fullPath)
        singleSegInps = []
        doubleSegInps = []

        # Sort it so that it completes the individual segments first
        for inp in validInpsList:
            if (inp.count('seg') == 1):
                singleSegInps.append(inp)
        for inp in validInpsList:
            if (inp.count('seg') == 2):
                doubleSegInps.append(inp)

        validInps = sorted(singleSegInps)

        monomersInSegment = []

        for orcaInput in validInps:
            print "Determining length of the segment:", orcaInput
            segmentNumbers = obtainSegmentNumbers(orcaInput[:-4])
            lengthOfSegment = determineMonomersInSegment(fullPath+orcaInput)
            monomersInSegment.append([segmentNumbers[0], lengthOfSegment])

        csvFileName = fullPath+'SegmentLengths.csv'

        csvFile = open(csvFileName, 'w+')
        csvWriter = csv.writer(csvFile, delimiter = ',')
        for length in monomersInSegment:
            csvWriter.writerow([length[0], length[1]])
        csvFile.close()

        print "Lengths of each segment written to", csvFileName
            
