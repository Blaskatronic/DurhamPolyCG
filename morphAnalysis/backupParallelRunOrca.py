import os
import numpy as np
import pylab as P
import time as T
import csv
import subprocess
import copy
import argparse

def getInpList(direc):
    fileList = os.listdir(direc)
    inpFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.inp'):
            inpFiles.append(str(fileName))
    return inpFiles

def sgesPresent(direc):
    fileList = os.listdir(direc)
    for fileName in fileList:
        if ".sge." in fileName:
            return True
    return False

def configureOrcaSubmission(orcaInput, fullPath, orcaPath, outputFiles, queue, parallel=True):
    templateHandle = open('./templates/serial.sge')
    submissionTemplate = templateHandle.readlines()
    templateHandle.close()

    outputPath, outputFiles = getOutputPath(fullPath, outputFiles)
    if outputPath == 0:
        return 1, outputFiles

    if parallel == False:
        runOrcaSerial = str(orcaPath)+'/orca '+str(fullPath)+' > '+str(outputPath)[:-4]+'.out'
        # print "Executing Orca..."
        subprocess.call(runOrcaSerial, shell=True)
        return 0, outputFiles

    if parallel == True:
        # print "Submitting ORCA..."
        submissionTemplate[-1] = str(orcaPath)+'/orca '+str(fullPath)+' > '+str(outputPath)[:-4]+'.out\n'
        SGEName = './sges/'+str(orcaInput)[:-4]+'.sge'
        newRun = open(SGEName, 'w+')

        # print "Writing SGE file as:", str(SGEName)
        newRun.writelines(submissionTemplate)
        newRun.close()
        if queue == "seq6.q":
            submissionOutput = os.popen("qsub -q "+queue+" "+SGEName[2:])
            submissionOutput = submissionOutput.read()
            print str(submissionOutput[:-1]), "for", fullPath
            jobNumber = int(str(submissionOutput).split(' ')[2])
        elif queue == "par6.q":
            submissionOutput = os.popen("qsub -q "+queue+" -pe orte 1 "+SGEName[2:])
            submissionOutput = submissionOutput.read()
            print str(submissionOutput[:-1]), "for", fullPath
            jobNumber = int(str(submissionOutput).split(' ')[2])
        else:
            print "Queue =", queue
            raise SystemError('Queue not correctly hardcoded')
        return jobNumber, outputFiles


def getOutputPath(inputPath, outputFiles=None):
    slashes = findIndex(inputPath, '/')
    test = inputPath.split('/')
    test[1] = 'orcaOutputs'
    outputPath = ''
    for element in test:
        outputPath += str(element)
        outputPath += '/'
    outputFileName = test[-1][:-3]+'out'
    if outputFiles == None:
        morphologyName = outputPath[slashes[1]+2:slashes[2]+1]
        outputMorphologies = os.listdir('./orcaOutputs/')
        if morphologyName not in outputMorphologies:
            os.makedirs('./orcaOutputs/'+morphologyName)
        outputFiles = os.listdir('./orcaOutputs/'+morphologyName)
    if outputFileName in outputFiles:
        print "\rFile", outputFileName, "already exists, continuing...",
        return 0, outputFiles
    return outputPath[:-1], outputFiles



def checkRunningJobs():
    runningJobs = 0
    queuedWaitingJobs = 0
    qstatOutput = subprocess.Popen("qstat", stdout=subprocess.PIPE).communicate()[0]
    qstatOutput = qstatOutput.split('\n')
    for lineNo in range(len(qstatOutput)):
        if "seg" in qstatOutput[lineNo]:
            # Make this more intelligent so that it can distinguish between running, queued etc.
            jobLine = []
            tempLine = qstatOutput[lineNo].split(' ')
            for element in tempLine:
                if (len(element) != 0) and (element != "\n"):
                    jobLine.append(element)
            if jobLine[4] == 'r':
                runningJobs += 1
            elif jobLine[4] == 'qw':
                queuedWaitingJobs += 1
            elif jobLine[4] == 'Eqw':
                print "CAUTION: ONE JOB HAS ERRORED OUT! PLEASE TERMINATE AND CHECK"
                raise SystemError('Job Submission Error')
    return runningJobs, queuedWaitingJobs




def findIndex(string, target):
    '''This function returns the locations of an inputted character (target) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == target:
            locations.append(index)
        index += 1
    return locations

def splitList(listToSplit, splitNumber):
    '''This function splits a list (listToSplit) into equal chunks of size splitNumber (and then yields any remaining elements)'''
    newList = []
    for i in range(0, len(listToSplit), splitNumber):
        newList.append(listToSplit[i:i+splitNumber])
    return newList


if __name__ == "__main__":

    # validMorphs = os.listdir('./orcaInputs')
    # validMorphs = sorted(validMorphs)

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--morphology", help="The morphology directory in ./orcaOutputs within which to analyse the ORCA .out files")
    parser.add_argument("-c", "--cores", help="The total number of simultaneous ORCA ZINDO calculations to perfrom on Hamilton")
    parser.add_argument("-q", "--queue", help="The queue to submit the hamilton jobs to")
    args = parser.parse_args()

    orcaPath = '/gpfs/scratch/ghsk28/ORCA/orca_3_0_3_linux_x86-64'

    completedInputFiles = 0
    totalInputFiles = 0
    t0 = T.time()

    morphology = args.morphology
    coresToUse = int(args.cores)
    queue = args.queue

    # for morphology in validMorphs:
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

    if len(validInps) == 0:
        for inp in validInpsList:
            singleSegInps.append(inp)
        validInps = sorted(singleSegInps)

    if totalInputFiles == 0:
        totalInputFiles = len(validInps)

    print "There are", totalInputFiles, "ORCA input files to run."


    t3 = T.time()

    totalSubmittedJobs = 0
    outputFiles = None

    if coresToUse > 0:
        while True:
            # Check how many jobs running
            if totalSubmittedJobs == len(validInps):
                print "All jobs submitted to Hamilton."
                break
            runningJobs, queuedJobs = checkRunningJobs()
            if (runningJobs + queuedJobs) < coresToUse:
                errorCode, outputFiles = configureOrcaSubmission(validInps[totalSubmittedJobs], fullPath+validInps[totalSubmittedJobs], orcaPath, outputFiles, queue, parallel=True)
                totalSubmittedJobs += 1
                # Clean up the blank output files that Hamilton generates
                removeSGEs = sgesPresent('./')
                if removeSGEs == True:
                    os.system("rm *.sge.*")
            else:
                T.sleep(2)

    else:
        for orcaInput in validInps:
            t1 = T.time()
            errorCode, outputFiles = configureOrcaSubmission(orcaInput, fullPath+validInps[completedInputFiles], orcaPath, outputFiles, queue, parallel=False)
            t2 = T.time()
            completedInputFiles += 1

            currentElapsedTime = float(t2)-float(t0)
            if currentElapsedTime < 60:
                timeUnits2 = 'seconds.'
            elif currentElapsedTime < 3600:
                seconds = str(currentElapsedTime%60)
                while len(seconds) > 4:
                    seconds = seconds[:-1]
                currentElapsedTime /= 60.0
                timeUnits2 = str('minutes, '+str(seconds)+' seconds.')
            elif currentElapsedTime < 86400:
                minutes = currentElapsedTime%60
                minutesString = str(int(minutes))
                seconds = str(minutes%60)
                while len(seconds) > 4:
                    seconds = seconds[:-1]
                currentElapsedTime /= 3600.0
                timeUnits2 = str('hours, '+str(minutesString)+' minutes and '+str(seconds)+' seconds.')
            averageTimePerInp = (float(t2)-float(t0))/float(completedInputFiles)
            expectedTotalTime = averageTimePerInp*totalInputFiles
            remainingTime = expectedTotalTime - currentElapsedTime
            if remainingTime < 60:
                timeunits = 'seconds.'
            elif remainingTime < 3600:
                seconds = str(currentElapsedTime%60)
                while len(seconds) > 4:
                    seconds = seconds[:-1]
                remainingTime /= 60.0
                timeunits = str('minutes, '+str(seconds)+' seconds.')
            elif remainingTime < 86400:
                minutes = currentElapsedTime%60
                minutesString = str(int(minutes))
                seconds = str(minutes%60)
                while len(seconds) > 4:
                    seconds = seconds[:-1]
                remainingTime /= 3600.0
                timeunits = str('hours, '+str(minutesString)+' minutes and '+str(seconds)+' seconds.')
            else:
                remainingTime /= 86400.0
                timeunits = 'days.'
            if errorCode == 0:
                print "\rORCA run completed for", orcaInput, "in %.1f seconds. Elasped time: %.1f %s Remaining time: %.1f %s" % (float(t2-t1), float(currentElapsedTime), str(timeUnits2), float(remainingTime), str(timeunits)),


    t4 = T.time()

    print "Calculations complete in %.3f seconds." % (float(t4)-float(t3))
