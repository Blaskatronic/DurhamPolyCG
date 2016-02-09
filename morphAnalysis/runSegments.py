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


def getFilesList(direc):
    fileList = os.listdir(direc)
#    files = []
#    for element in fileList:
#        files.append(element[7:])
    datFiles = []
    lammpsFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.dat'):
            datFiles.append(str(fileName))
    popList = []
    for datNumber in range(len(datFiles)):
        for fileName in fileList:
            addedLammps = 0
            if (datFiles[datNumber][:-4] in fileName) and ('.lammpstrj' in fileName):
                lammpsFiles.append(fileName)
                addedLammps = 1
                break
        if addedLammps == 0:
            popList.append(datNumber)
    popList = sorted(popList, reverse=True)
    for popIndex in popList:
        datFiles.pop(popIndex)
    
    return datFiles, lammpsFiles


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
   

def parallelSort(list1, list2):
    data = zip(list1, list2)
    data.sort()
    list1, list2 = map(lambda t: list(t), zip(*data))
    return list1, list2


if __name__ == '__main__':
    datFiles, lammpstrjFiles = getFilesList('./')
    datFiles, lammpstrjFiles = parallelSort(datFiles, lammpstrjFiles)
    usableFiles = zip(datFiles, lammpstrjFiles)
    exitFlag = 0
    while exitFlag == 0:
        while True:
            print "\n---=== FILES THAT CAN BE RUN BY SEGMENT.PY ===---"
            for elementNo in range(len(usableFiles)):
                print str(elementNo)+"):", usableFiles[elementNo]
            print str(elementNo+1)+"): Exit runSegments.py"
            # print "Valid files =", zip(datFiles, lammpstrjFiles)
            runThisFile = raw_input("Please pick a file to run (integer, default = 0): ")
            if len(runThisFile) == 0:
                runThisFile = 0
            else:
                try:
                    runThisFile = int(runThisFile)
                except:
                    print "Please enter an integer between 0 and", len(usableFiles)
                    continue
            if (runThisFile < 0) or (runThisFile > len(usableFiles)):
                print "Please enter an integer between 0 and", len(usableFiles)
                continue
            elif runThisFile == len(usableFiles):
                print "Exiting Program..."
                exitFlag = 1
                break
            break
        if exitFlag == 0:
            print "Checking directories are in place for", str(usableFiles[runThisFile])+"..."
            checkDirectories(usableFiles[runThisFile][0])
            t0 = T.time()
            subprocess.call('python segments.py -d '+str(usableFiles[runThisFile][0])+' -l '+str(usableFiles[runThisFile][1]), shell=True)
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
            print "Segments.py calculations complete in %.1f %s." % (float(elapsedTime), str(timeunits))
            print "----------====================----------"

    print "Battle control terminated."
