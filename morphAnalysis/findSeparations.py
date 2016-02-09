import os
import numpy as np
import pylab as P
import random as R
import csv
import mpl_toolkits.mplot3d.axes3d as p3
import time as T
import copy
import csv
import argparse
import math
import subprocess


def getEquilibratedMorphologies(direc):
    fileList = os.listdir(direc)
#    files = []
#    for element in fileList:
#        files.append(element[7:])
    equilibratedMorphologies = []
    for fileName in fileList:
        if ("Equil" in fileName):
            equilibratedMorphologies.append(str(fileName))
    return equilibratedMorphologies


def checkDirectories(datFileName):
    dotsLoc = findIndex(datFileName, '.')
    actualName = datFileName[:dotsLoc[-1]]
    inputDir = os.listdir('./orcaInputs')
    outputDir = os.listdir('./orcaOutputs')
    if (actualName not in inputDir):
        if (actualName not in outputDir):
            print "Neither input nor output directories present!"
            os.makedirs('./orcaInputs/'+str(actualName))
            os.makedirs('./orcaOutputs/'+str(actualName))
            print "Correct directories created."
        else:
            print "Input directory not present!"
            os.makedirs('./orcaInputs/'+str(actualName))
            print "Correct directory created."
    else:
        print "Both input and output directories present and correct."


def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    return locations
   

def readSingleSegmentsCSV(morphology):
    singleSegments = {}
    CSVpath = './orcaOutputs/'+str(morphology)+'/SingleSegments.csv'
    CSVFile = open(CSVpath, 'r')
    CSVreader = csv.reader(CSVFile, delimiter=',')
    for row in CSVreader:
        segmentID = row[0]
        segmentCoords = [row[3], row[4], row[5]]
        singleSegments[segmentID] = segmentCoords
    CSVFile.close()
    return singleSegments


def readTransferIntegralsCSV(morphology):
    dualSegments = []
    CSVpath = './orcaOutputs/'+str(morphology)+'/TransferIntegrals.csv'
    CSVFile = open(CSVpath, 'r')
    CSVreader = csv.reader(CSVFile, delimiter=',')
    for row in CSVreader:
        segment1ID = row[0]
        segment2ID = row[1]
        dualSegments.append([segment1ID, segment2ID])
    CSVFile.close()
    return dualSegments
    



def findSeparations(morphology):
    singleSegments = readSingleSegmentsCSV(morphology)
    dualSegments = readTransferIntegralsCSV(morphology)
    # GET THE SIM VOL DATA


def parallelSort(list1, list2):
    data = zip(list1, list2)
    data.sort()
    list1, list2 = map(lambda t: list(t), zip(*data))
    return list1, list2


if __name__ == '__main__':
    morphologies = getEquilibratedMorphologies('./orcaOutputs')
    exitFlag = 0
    while exitFlag == 0:
        while True:
            print "\n---=== MORPHOLOGIES THAT WE CAN CALCULATE SEPARATION DATA FOR ===---"
            print "IMPORTANT NOTE: THIS PROGRAM DOES NOT TAKE INTO ACCOUNT PERIODIC HOPS SO THE SEPARATIONS FOR THESE PAIRS WILL BE MUCH LARGER THAN EXPECTED"
            if len(morphologies) == 0:
                print "ERROR: No ORCA Outputs found in ./orcaOutputs. Please runOrca.py first!"
                exitFlag = 1
                break
            for elementNo in range(len(morphologies)):
                print str(elementNo)+"):", morphologies[elementNo]
            print str(elementNo+1)+"): Exit findSeparations.py"
            # print "Valid files =", zip(datFiles, lammpstrjFiles)
            runThisFile = raw_input("Please pick a file to run (integer, default = 0): ")
            if len(runThisFile) == 0:
                runThisFile = 0
            else:
                try:
                    runThisFile = int(runThisFile)
                except:
                    print "Please enter an integer between 0 and", len(morphologies)
                    continue
            if (runThisFile < 0) or (runThisFile > len(morphologies)):
                print "Please enter an integer between 0 and", len(morphologies)
                continue
            elif runThisFile == len(morphologies):
                print "Exiting Program..."
                exitFlag = 1
                break
            break
        if exitFlag == 0:
            print "\n"
            t0 = T.time()
            findSeparations(morphologies[runThisFile])
            subprocess.call('python analyseOrca.py -m '+str(morphologies[runThisFile]), shell=True)
            t1 = T.time()
            elapsedTime = float(t1) - float(t0)
            if elapsedTime < 60:
                timeunits = 'seconds.'
            elif elapsedTime < 3600:
                elapsedTime /= 60.0
                timeunits = 'minutes.'
            elif elapsedTime < 86400:
                elapsedTime /= 3600.0
                timeunits = 'hours.'
            else:
                elapsedTime /= 86400.0
                timeunits = 'days.'
            print "----------====================----------"
            print "findSeparations.py calculations complete in %.1f %s." % (float(elapsedTime), str(timeunits))
            print "----------====================----------"

    print "Battle control terminated."
