import os
import numpy as np
import pylab as P
import time as T
import csv
import subprocess
import copy
import argparse
import sys
import multiprocessing as mp
import itertools

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
        if ".sge.o" in fileName:
            return True
    return False



def testRun(orcaInput, fullPath, orcaPath, outputPath):
    runOrcaCommand = str(orcaPath)+'/orca '+str(fullPath)+' > '+str(outputPath)[:-4]+'.out'
    print runOrcaCommand
    subprocess.call(runOrcaCommand, shell=True)



def configureOrcaSubmission(orcaInput, fullPath, orcaPath, outputFiles, engine, queue, parallel=True):
    outputPath, outputFiles = getOutputPath(fullPath, outputFiles)
    if outputPath == 0:
        return 0, outputFiles, None

    if parallel == False:
        runOrcaSerial = str(orcaPath)+'/orca '+str(fullPath)+' > '+str(outputPath)[:-4]+'.out'
        # print "Executing Orca..."
        subprocess.call(runOrcaSerial, shell=True)
        return 0, outputFiles, None

    if parallel == True:
        if engine == 'sge':
            templateHandle = open('./templates/serial.sge')
        elif engine == 'slurm':
            templateHandle = open('./templates/serial.slurm')
        submissionTemplate = templateHandle.readlines()
        templateHandle.close()
        # print "Submitting ORCA..."
        submissionTemplate[-1] = str(orcaPath)+'/orca '+str(fullPath)+' > '+str(outputPath)[:-4]+'.out\n'

        if engine == 'sge':
            SGEName = './'+str(orcaInput)[:-4]+'.sge'
        elif engine == 'slurm':
            SGEName = './'+str(orcaInput)[:-4]+'.slurm'
        newRun = open(SGEName, 'w+')

        # print "Writing SGE file as:", str(SGEName)
        newRun.writelines(submissionTemplate)
        newRun.close()
        if (queue == "seq6.q") or (queue == "seq5.q"):
            while True:
                if engine == 'sge':
                    submissionOutput = os.popen("qsub -q "+queue+" "+SGEName[2:])
                elif engine == 'slurm':
                    submissionOutput = os.popen("sbatch -p "+queue+" -n 1 "+SGEName[2:])
                submissionOutput = submissionOutput.read()
                if ("failed receiving gdi request" in submissionOutput) or ("sbatch: error:" in submissionOutput):
                    continue
                else:
                    break
            print str(submissionOutput[:-1]), "to", str(engine)+"/"+str(queue), "for", fullPath
            if engine == 'sge':
                jobNumber = int(str(submissionOutput).split(' ')[2])
            elif engine == 'slurm':
                jobNumber = int(str(submissionOutput).split(' ')[3])
        elif (queue == "par6.q") or (queue == "par5.q") or (queue == "par4.q"):
            while True:
                if engine == 'sge':
                    submissionOutput = os.popen("qsub -q "+queue+" -pe orte 1 "+SGEName[2:])
                    submissionOutput = submissionOutput.read()
                    jobNumber = int(str(submissionOutput).split(' ')[2])
                elif engine == 'slurm':
                    submissionOutput = os.popen("sbatch -p "+queue+" -n 1 --ntasks-per-node=1 --cpus-per-task=1 "+SGEName[2:])
                    submissionOutput = submissionOutput.read()
                    jobNumber = int(str(submissionOutput).split(' ')[3])
                if ("failed receiving gdi request" in submissionOutput) or ("sbatch: error:" in submissionOutput):
                    continue
                else:
                    break
            print str(submissionOutput[:-1]), "to", str(engine)+"/"+str(queue), "for", fullPath
        else:
            print "Queue =", queue
            raise SystemError('Queue not correctly hardcoded')
        os.system("rm "+str(SGEName[2:]))
        return jobNumber, outputFiles, SGEName


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



def checkRunningJobs(engine, queue):
    runningJobs = 0
    queuedWaitingJobs = 0
    if engine == 'sge':
        runningOutput = subprocess.Popen("qstat", stdout=subprocess.PIPE).communicate()[0]
        runningOutput = runningOutput.split('\n')
    elif engine == 'slurm':
        squeueCommand = ['squeue', '-p', queue]
        runningOutput = subprocess.Popen(squeueCommand, stdout=subprocess.PIPE).communicate()[0]
        runningOutput = runningOutput.split('\n')
    runningJobNumbers = []
    for lineNo in range(len(runningOutput)):
        if "seg" in runningOutput[lineNo]:
            # Make this more intelligent so that it can distinguish between running, queued etc.
            jobLine = []
            tempLine = runningOutput[lineNo].split(' ')
            for element in tempLine:
                if (len(element) != 0) and (element != "\n"):
                    jobLine.append(element)
            if engine == 'sge':
                if jobLine[4] != 't':
                    runningJobNumbers.append(int(jobLine[0]))
                if jobLine[4] == 'r':
                    runningJobs += 1
                elif jobLine[4] == 'qw':
                    queuedWaitingJobs += 1
                elif jobLine[4] == 'Eqw':
                    print "CAUTION: ONE JOB HAS ERRORED OUT! PLEASE TERMINATE AND CHECK"
                    raise SystemError('Job Submission Error')
            elif engine == 'slurm':
                if jobLine[4] != 'T':
                    runningJobNumbers.append(int(jobLine[0]))
                if jobLine[4] == 'R':
                    runningJobs += 1
                elif jobLine[4] == 'PD':
                    queuedWaitingJobs += 1
    return runningJobs, queuedWaitingJobs, runningJobNumbers




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
    parser.add_argument("-e", "--engine", help="The Hamilton submission engine to use for job submission")
    parser.add_argument("-q", "--queue", help="The queue to submit the hamilton jobs to")
    args = parser.parse_args()

    orcaPath = '/gpfs/scratch/ghsk28/ORCA/orca_3_0_3_linux_x86-64'

    completedInputFiles = 0
    totalInputFiles = 0
    t0 = T.time()

    morphology = args.morphology
    coresToUse = int(args.cores)
    queue = args.queue
    engine = args.engine

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

    jobsRunning = 0
    jobsToWaitFor = []
    completedInputFiles = 0
    outputFiles = None


    # output

    # procPool = mp.Pool(processes=coresToUse)
    # procPool.map(runStart, itertools.izip(validInps, itertools.repeat(fullPath), itertools.repeat(outputPath))


    coresToUse = mp.cpu_count()

    print "CORES TO USE =", coresToUse



    for orcaInput in validInps:
        printWaiting = False
        printNewLine = False

        outputPath, outputFiles = getOutputPath(fullPath+orcaInput, outputFiles)
        if outputPath == 0:
            continue

        orcaProcess = mp.Process(target=testRun, args=(orcaInput, fullPath+orcaInput, orcaPath, outputPath))
        orcaProcess.start()
        jobsRunning += 1
        jobsToWaitFor.append(orcaProcess)
        if jobsRunning >= coresToUse:
            popList = []
            jobFinished = False
            while jobFinished == False:
                for jobNo in range(len(jobsToWaitFor)):
                    if jobsToWaitFor[jobNo].is_alive() == False:
                        popList.append(jobNo)
                        completedInputFiles += 1
                        jobsRunning -= 1
                        jobFinished = True
                popList = sorted(popList, reverse=True)
                for popIndex in popList:
                    jobsToWaitFor.pop(popIndex)
                if jobFinished == False:
                    T.sleep(2)

    # raise SystemError('CUSTARD')





















    # for orcaInput in validInps:
    #     t1 = T.time()
    #     if coresToUse == 0:
    #         parallel = False
    #         jobNumber, outputFiles, SGEName = configureOrcaSubmission(orcaInput, fullPath+orcaInput, orcaPath, outputFiles, engine, queue, parallel)
    #         completedInputFiles += 1
    #         t2 = T.time()

    #         currentElapsedTime = float(t2)-float(t0)
    #         if currentElapsedTime < 60:
    #             timeUnits2 = 'seconds.'
    #         elif currentElapsedTime < 3600:
    #             seconds = str(currentElapsedTime%60)
    #             while len(seconds) > 4:
    #                 seconds = seconds[:-1]
    #             currentElapsedTime /= 60.0
    #             timeUnits2 = str('minutes, '+str(seconds)+' seconds.')
    #         elif currentElapsedTime < 86400:
    #             minutes = currentElapsedTime%60
    #             minutesString = str(int(minutes))
    #             seconds = str(minutes%60)
    #             while len(seconds) > 4:
    #                 seconds = seconds[:-1]
    #             currentElapsedTime /= 3600.0
    #             timeUnits2 = str('hours, '+str(minutesString)+' minutes and '+str(seconds)+' seconds.')
    #         averageTimePerInp = (float(t2)-float(t0))/float(completedInputFiles)
    #         expectedTotalTime = averageTimePerInp*totalInputFiles
    #         remainingTime = expectedTotalTime - currentElapsedTime
    #         if remainingTime < 60:
    #             timeunits = 'seconds.'
    #         elif remainingTime < 3600:
    #             seconds = str(currentElapsedTime%60)
    #             while len(seconds) > 4:
    #                 seconds = seconds[:-1]
    #             remainingTime /= 60.0
    #             timeunits = str('minutes, '+str(seconds)+' seconds.')
    #         elif remainingTime < 86400:
    #             minutes = currentElapsedTime%60
    #             minutesString = str(int(minutes))
    #             seconds = str(minutes%60)
    #             while len(seconds) > 4:
    #                 seconds = seconds[:-1]
    #             remainingTime /= 3600.0
    #             timeunits = str('hours, '+str(minutesString)+' minutes and '+str(seconds)+' seconds.')
    #         else:
    #             remainingTime /= 86400.0
    #             timeunits = 'days.'
    #         print "\rORCA run completed for", orcaInput, "in %.1f seconds. Elasped time: %.1f %s Remaining time: %.1f %s" % (float(t2-t1), float(currentElapsedTime), str(timeUnits2), float(remainingTime), str(timeunits)),

    #     else:
    #         printWaiting = False
    #         printNewLine = False
    #         testJobs = []
    #         while True:
    #                 jobNumber, outputFiles, SGEName = configureOrcaSubmission(orcaInput, fullPath+orcaInput, orcaPath, outputFiles, engine, queue, parallel)
    #                 if jobNumber == 0:
    #                     # File already exists, so skip this one
    #                     printNewLine = True
    #                     break
    #                 else:
    #                     if printNewLine == True:
    #                         printNewLine = False
    #                         print "\n"
    #                 jobsRunning += 1
    #                 jobsToWaitFor.append(jobNumber)
    #                 printWaiting = True
    #                 if engine == 'sge':
    #                     removeSGEs = sgesPresent('./')
    #                     if removeSGEs == True:
    #                         os.system("rm *.sge.o*")
    #                 break
    #             # else:
    #             #     if printWaiting == True:
    #             #         print "Waiting for submission slot..."
    #             #     printWaiting = False
    #             #     T.sleep(2)
    #             #     runningJobs, queuedJobs, runningJobNumbers = checkRunningJobs(engine, queue)
    #             #     for jobNumber in jobsToWaitFor:
    #             #         if jobNumber not in runningJobNumbers:
    #             #             jobsToWaitFor.remove(jobNumber)
    #             #             jobsRunning -= 1
    #             #             completedInputFiles += 1
    #             #             t2 = T.time()


    # t4 = T.time()
    # print "Calculations complete in %.3f seconds." % (float(t4)-float(t3))


