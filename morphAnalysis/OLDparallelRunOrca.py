import os
import numpy as np
import pylab as P
import time as T
import csv
import subprocess
import copy
import multiprocessing

def getInpList(direc):
    fileList = os.listdir(direc)
    inpFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.inp'):
            inpFiles.append(str(fileName))
    return inpFiles

def configureOrcaSubmission(orcaInput, fullPath, orcaPath):
    templateHandle = open('./templates/serial.sge')
    submissionTemplate = templateHandle.readlines()
    templateHandle.close()

    outputPath = getOutputPath(fullPath)
    if outputPath == 0:
        return 1

    submissionTemplate[-1] = str(orcaPath)+'/orca '+str(fullPath)+' > '+str(outputPath)[:-4]+'.out\n'
    SGEName = './'+str(orcaInput)[:-4]+'.sge'
    newRun = open(SGEName, 'w+')

    # print "Writing SGE file as:", str(SGEName)
    # newRun.writelines(submissionTemplate)
    # newRun.close()

    runOrca = submissionTemplate[-1][:-1]

    print "Executing Orca..."
    subprocess.call(runOrca, shell=True)
    return 0


    # if procSlots > 1:
    #     os.system("qsub -q par6.q -pe orte "+str(procSlots)+" "+SGEName[2:])
    # else:
    #     os.system("qsub -q seq6.q "+SGEName[2:])


def getOutputPath(inputPath):
    slashes = findIndex(inputPath, '/')
    test = inputPath.split('/')
    test[1] = 'orcaOutputs'
    outputPath = ''
    for element in test:
        outputPath += str(element)
        outputPath += '/'
    outputFileName = test[-1][:-3]+'out'
    morphologyName = outputPath[slashes[1]+2:slashes[2]+1]
    outputMorphologies = os.listdir('./orcaOutputs/')
    if morphologyName not in outputMorphologies:
        os.makedirs('./orcaOutputs/'+morphologyName)
    outputFiles = os.listdir('./orcaOutputs/'+morphologyName)
    if outputFileName in outputFiles:
        print "File", outputFileName, "already exists, continuing...\n"
        return 0
    return outputPath[:-1]


def findIndex(string, target):
    '''This function returns the locations of an inputted character (target) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == target:
            locations.append(index)
        index += 1
    return locations


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

        validInps = sorted(singleSegInps)+sorted(doubleSegInps)

        if totalInputFiles == 0:
            totalInputFiles = len(validMorphs)*len(validInps)

        parallelProcesses = []


        t3 = T.time()

        for orcaInput in validInps:
            print "Running ORCA for file:", orcaInput
            t1 = T.time()
            orcaProcess = multiprocessing.Process(target = configureOrcaSubmission, args = (orcaInput, fullPath+orcaInput, '/gpfs/scratch/ghsk28/ORCA/orca_3_0_3_linux_x86-64'), name = str(orcaInput))
            # errorCode = configureOrcaSubmission(orcaInput, fullPath+orcaInput, '/gpfs/scratch/ghsk28/ORCA/orca_3_0_3_linux_x86-64')
            orcaProcess.start()
            parallelProcesses.append(orcaProcess)
            t2 = T.time()
            completedInputFiles += 1

            currentElapsedTime = t2-t0
            averageTimePerInp = (t2-t0)/float(completedInputFiles)
            expectedTotalTime = averageTimePerInp*totalInputFiles
            remainingTime = expectedTotalTime - currentElapsedTime
            if remainingTime < 60:
                timeunits = 'seconds.'
            elif remainingTime < 3600:
                remainingTime /= 60.0
                timeunits = 'minutes.'
            elif remainingTime < 86400:
                remainingTime /= 3600.0
                timeunits = 'hours.'
            else:
                remainingTime /= 86400.0
                timeunits = 'days.'


            # if errorCode == 0:
            #     print "ORCA run completed in %.1f seconds. Remaining time: %.1f %s" % (float(t2-t1), float(remainingTime), str(timeunits))


        for orcaProcess in parallelProcesses:
            print "ZINDO calculations for", orcaProcess.name, "are now complete."
            orcaProcess.join()

        t4 = T.time()

        print "Calculations complete in %.3f seconds." % (float(t4)-float(t3))

        print "Finished?"

