import os
import sys
import numpy as np
import pylab as P
import time as T
import csv
import subprocess
import argparse
import copy
import multiprocessing as mp

def getOutList(direc):
    fileList = os.listdir(direc)
    inpFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.out'):
            inpFiles.append(str(fileName))
    return inpFiles

def getOrbitals(outputFilePath):
    ORCAfile = open(outputFilePath, 'r')
    ORCALines = ORCAfile.readlines()
    ORCAfile.close()
    molecularOrbitals = []

    recordingLines = False
    startRecordingOnLine = 1E20
    for lineNo in range(len(ORCALines)):
        if lineNo == startRecordingOnLine:
            recordingLines = True
        if ("SCF NOT CONVERGED" in ORCALines[lineNo]) or ("SERIOUS PROBLEM" in ORCALines[lineNo]):
            print "WARNING: SCF was not converged for file", str(outputFilePath)+"."
            print "No molecular orbital data exists."
            print "Rerun/modify morphology to obtain molecular orbital data"
            return False
        if "ORBITAL ENERGIES" in ORCALines[lineNo]:
            startRecordingOnLine = lineNo+4
        if (recordingLines == True) and ("------" in ORCALines[lineNo]):
            recordingLines = False
            break
        # Actually record data
        if (recordingLines == True):
            splitLine = ORCALines[lineNo].split(' ')
            molecularOrbital = []
            firstElement = True
            for element in splitLine:
                if (len(element) != 0) and (element != '\n'):
                    if firstElement == True:
                        molecularOrbital.append(int(element))
                        firstElement = False
                    else:
                        molecularOrbital.append(float(element))
            molecularOrbitals.append(molecularOrbital)
    return molecularOrbitals


def getAtomPositions(outputFilePath):
    ORCAfile = open(outputFilePath, 'r')
    ORCALines = ORCAfile.readlines()
    ORCAfile.close()
    atomPositions = []

    recordingLines = False
    startRecordingOnLine = 1E20
    for lineNo in range(len(ORCALines)):
        if lineNo == startRecordingOnLine:
            recordingLines = True
        if "CARTESIAN COORDINATES (ANGSTROEM)" in ORCALines[lineNo]:
            startRecordingOnLine = lineNo+2
        if (recordingLines == True) and (ORCALines[lineNo] == '\n'):
            recordingLines = False
            break
        # Actually record data
        if (recordingLines == True):
            splitLine = ORCALines[lineNo].split(' ')
            atom = []
            firstElement = True
            for element in splitLine:
                if (len(element) != 0) and (element != '\n'):
                    if firstElement == True:
                        atom.append(str(element))
                        firstElement = False
                    else:
                        atom.append(float(element))
            atomPositions.append(atom)
    return atomPositions


def determineCOM(inputAtoms):
    massWeightedX = 0.
    massWeightedY = 0.
    massWeightedZ = 0.
    totalMass = 0.
    for atom in inputAtoms:
        if atom[0] == 0: # Determining chain COM, therefore mass factor is irrelevent
            mass = 1.0                
        elif atom[0] == 'S':
            mass = 32.065
        elif atom[0] == 'C':
            mass = 12.0107
        elif atom[0] == 'H':
            mass = 1.00794

        massWeightedX += atom[1]*mass
        massWeightedY += atom[2]*mass
        massWeightedZ += atom[3]*mass
        totalMass += mass
    return np.array([massWeightedX/float(totalMass), massWeightedY/float(totalMass), massWeightedZ/float(totalMass)])


def obtainSplittingEnergies(orbitals):
    for orbitalNo in range(len(orbitals)):
        if orbitals[orbitalNo][1] == 0.0:
            LUMONo = orbitalNo
            break
    HOMO = orbitals[LUMONo-1][3] # The HOMO level in eV
    HOMO_1 = orbitals[LUMONo-2][3] # The HOMO-1 level in eV
    return [HOMO, HOMO_1]
        
        


def obtainImportantEnergies(segmentNumbers, orbitals):
    HOMO, HOMO_1 = obtainSplittingEnergies(orbitals)
    if len(segmentNumbers) == 1:
        importantEnergies = [segmentNumbers[0], HOMO]
    elif len(segmentNumbers) == 2:
        importantEnergies = [segmentNumbers[0], segmentNumbers[1], HOMO, HOMO_1]
    else:
        print "Name =", name
        print "Number of Segments =", numberOfSegments
        raise SystemError('Number of Segments == 0 or > 2')  
    return importantEnergies


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


def turnOffSOSCF(outputFile, fullInputPath):
    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
    originalLines = fileName.readlines()
    fileName.close()

    originalLines[3] = '!ZINDO/S NoSOSCF\n'

    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
    fileName.writelines(originalLines)
    fileName.close()


def reduceTolerance(outputFile, fullInputPath):
    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
    originalLines = fileName.readlines()
    fileName.close()

    originalLines[3] = '!ZINDO/S NoSOSCF SloppySCF\n'

    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
    fileName.writelines(originalLines)
    fileName.close()


def increaseIterations(outputFile, fullInputPath):
    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
    originalLines = fileName.readlines()
    fileName.close()

    originalLines.append('\n%scf MaxIter 500 end')

    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
    fileName.writelines(originalLines)
    fileName.close()


def increaseGrid(outputFile, fullInputPath):
    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
    originalLines = fileName.readlines()
    fileName.close()

    originalLines[3] = '!ZINDO/S SlowConv Grid7 NoFinalGrid\n'

    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
    fileName.writelines(originalLines)
    fileName.close()


def increaseGridNoSOSCF(outputFile, fullInputPath):
    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
    originalLines = fileName.readlines()
    fileName.close()

    originalLines[3] = '!ZINDO/S SlowConv Grid7 NoFinalGrid NoSOSCF SloppySCF\n'
    originalLines.append('\n%scf MaxIter 500 end')

    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
    fileName.writelines(originalLines)
    fileName.close()


def revertORCAFiles(outputFile, fullInputPath):
    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
    originalLines = fileName.readlines()
    fileName.close()
    originalLines[3] = '! ZINDO/S\n'
    for lineNo in range(len(originalLines)):
        # REMOVE THE SCF ITER
        if "%scf MaxIter" in originalLines[lineNo]:
            originalLines.pop(lineNo)
            break

    fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
    fileName.writelines(originalLines)
    fileName.close()



def testRun(orcaInput, fullPath, orcaPath):
    outputPath = getOutputPath(fullPath)
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
                    submissionOutput = os.popen("sbatch -p "+queue+" "+SGEName[2:])
                submissionOutput = submissionOutput.read()
                if "failed receiving gdi request" in submissionOutput:
                    continue
                else:
                    break
            print str(submissionOutput[:-1]), "to", str(engine)+"/"+str(queue), "for", fullPath
            if engine == 'sge':
                jobNumber = int(str(submissionOutput).split(' ')[2])
            elif engine == 'slurm':
                jobNumber = int(str(submissionOutput).split(' ')[3])
        elif (queue == "par6.q") or (queue == "par5.q")or (queue == "par4.q"):
            if engine == 'sge':
                submissionOutput = os.popen("qsub -q "+queue+" -pe orte 1 "+SGEName[2:])
                submissionOutput = submissionOutput.read()
                jobNumber = int(str(submissionOutput).split(' ')[2])
            elif engine == 'slurm':
                submissionOutput = os.popen("sbatch -p "+queue+" -n 1 "+SGEName[2:])
                submissionOutput = submissionOutput.read()
                jobNumber = int(str(submissionOutput).split(' ')[3])
            print str(submissionOutput[:-1]), "to", str(engine)+"/"+str(queue), "for", fullPath
        else:
            print "Queue =", queue
            raise SystemError('Queue not correctly hardcoded')
        os.system("rm "+str(SGEName[2:]))
        return jobNumber, outputFiles, SGEName


def getOutputPath(inputPath, outputFiles=None):
    # print "InputPath =", inputPath
    slashes = findIndex(inputPath, '/')
    test = inputPath.split('/')
    test[1] = 'orcaOutputs'
    outputPath = ''
    for element in test:
        outputPath += str(element)
        outputPath += '/'
    outputFileName = test[-1][:-3]+'out'
    # print "outputFileName =", outputFileName
    # if outputFiles == None:
    #     morphologyName = outputPath[slashes[1]+2:slashes[2]+1]
    #     outputMorphologies = os.listdir('./orcaOutputs/')
    #     # print "outputMorphologies =", outputMorphologies
    #     # print "morphologyName =", morphologyName
    #     if morphologyName not in outputMorphologies:
    #         print "Making directory: ./orcaOutputs/"+morphologyName
    #         os.makedirs('./orcaOutputs/'+morphologyName)
    #     outputFiles = os.listdir('./orcaOutputs/'+morphologyName)

    # DON'T NEED THIS BIT - WE ALWAYS WANT TO RERUN IT AND OVERWRITE THE OLD ONE
    # if outputFileName in outputFiles:
    #     print "\rFile", outputFileName, "already exists, continuing...",
    #     return 0, outputFiles
    return outputPath[:-1]#, outputFiles




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

def sgesPresent(direc):
    fileList = os.listdir(direc)
    for fileName in fileList:
        if ".sge.o" in fileName:
            return True
    return False

def calculateTransferIntegral(seg1, seg2, E_HOMO, E_HOMO_1, epsilon1, epsilon2, koopmansApproximation=True):
    # print "E_HOMO - E_HOMO_1 =", E_HOMO-E_HOMO_1
    # print "epsilon1 - epsilon2 =", epsilon1-epsilon2
    if koopmansApproximation == False:
        if ((epsilon1 - epsilon2)**2 > (E_HOMO - E_HOMO_1)**2):
            return 0.5*(E_HOMO - E_HOMO_1), 1
        else:
            return 0.5*np.sqrt( (E_HOMO - E_HOMO_1)**2 - (epsilon1 - epsilon2)**2  ), 0
    else:
        return 0.5*(E_HOMO - E_HOMO_1), 0



# def calculateTransferIntegrals(HOMOEnergies, HOMOSplitting):
#     transferIntegrals = []
#     for pair in HOMOSplitting:
#         seg1 = pair[0]
#         seg2 = pair[1]
#         E_HOMO = pair[2]
#         E_HOMO_1 = pair[3]
#         epsilon1 = HOMOEnergies[seg1]
#         epsilon2 = HOMOEnergies[seg2]

#         transferIntegral = 0.5*np.sqrt( (E_HOMO - E_HOMO_1)**2 - (epsilon1 - epsilon2)**2  )

#         transferIntegrals.append([seg1, seg2, transferIntegral])
#     return transferIntegrals



        
def generateCSV(transferIntegrals, fullPath):
    fileName = fullPath+'TransferIntegrals.csv'
    # transferIntegralLines = []
    # csvReadOnly = csv.reader(open(fileName, 'r'), delimiter = ',')
    # for row in csvReadOnly:
    #     transferIntegralLines.append(row)


    csvFile = csv.writer(open(fileName, 'w+'), delimiter = ',')
    for transferIntegral in transferIntegrals:
        if transferIntegral not in transferIntegralLines:
            csvFile.writerow(transferIntegral)



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
    orcaPath = '/gpfs/scratch/ghsk28/ORCA/orca_3_0_3_linux_x86-64'

    failedFiles = []
    csvData = []


    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--morphology", help="The morphology directory in ./orcaOutputs within which to analyse the ORCA .out files")
    parser.add_argument("-c", "--cores", help="The total number of simultaneous ORCA ZINDO calculations to perfrom on Hamilton")
    parser.add_argument("-e", "--engine", help="The Hamilton submission engine to use for job submission")
    parser.add_argument("-q", "--queue", help="The queue to submit the hamilton jobs to")
    args = parser.parse_args()

    koopmansApproximation = True


    # validMorphs = os.listdir('./orcaOutputs')
    # validMorphs = sorted(validMorphs)

    HOMOEnergies = {}
    HOMOSplitting = {}
    segmentLengths = {}

    morphology = args.morphology
    coresToUse = int(args.cores)
    queue = args.queue
    engine = args.engine

    # for morphology in validMorphs:
    fullInputPath = './orcaInputs/'+morphology+'/'
    fullPath = './orcaOutputs/'+morphology+'/'
    assortedValidOuts = getOutList(fullPath)
    singleSegOuts = []
    doubleSegOuts = []
    for orcaOutput in assortedValidOuts:
        if orcaOutput.count("seg") == 1:
            singleSegOuts.append(orcaOutput)
        else:
            doubleSegOuts.append(orcaOutput)

    singleSegOuts = sorted(singleSegOuts)
    doubleSegOuts = sorted(doubleSegOuts)
    # validOuts = singleSegOuts + doubleSegOuts




    energyFileName = fullPath+'SingleSegments.csv'

    csvReadOnly = open(energyFileName, 'a+')
    csvReadOnly.seek(0)
    csvData = csv.reader(csvReadOnly, delimiter = ',')
    # if len(csvData) != 0:
    for row in csvData:
        HOMOEnergies[int(row[0])] = float(row[1])
    csvReadOnly.close()


    segmentLengthFileName = './orcaInputs/'+morphology+'/SegmentLengths.csv'

    csvReadOnly = open(segmentLengthFileName, 'r')
    csvData = csv.reader(csvReadOnly, delimiter = ',')
    for row in csvData:
        segmentLengths[int(row[0])] = int(row[1])
    csvReadOnly.close()


    # failedRerunCounter = 0
    # moreIter = False
    # previousLengthOfFailedSingles = 0
    failedSinglesDict = {} # A dictionary of the names of the failed singles and the number of failures for each
    outputFiles = None


    coresToUse = mp.cpu_count()
    

    print "CORES TO USE =", coresToUse




    while True:
        failedSingles = []
        for orcaOutput in singleSegOuts:
            segmentNumbers = obtainSegmentNumbers(orcaOutput[:-4])
            if (segmentNumbers[0] in HOMOEnergies) == False:
                print "Analysing ORCA single-segment file:", orcaOutput
                print "Segment not present in csv, analysing ORCA file..."
                molecularOrbitals = getOrbitals(fullPath+orcaOutput)
                if molecularOrbitals == False:
                    failedSingles.append(orcaOutput)
                    continue
                try:
                    importantEnergies = obtainImportantEnergies(segmentNumbers, molecularOrbitals)
                except:
                    failedSingles.append(orcaOutput)
                    continue
                # Remove the failure state for this output file if it exists
                failedSinglesDict.pop(orcaOutput, None)
                HOMOEnergies[segmentNumbers[0]] = importantEnergies[1]
                atomPositions = getAtomPositions(fullPath+orcaOutput)
                segmentCOM = determineCOM(atomPositions)
                csvFile = open(energyFileName, 'a')
                csvWriter = csv.writer(csvFile, delimiter = ',')
                # Write the CSV as [SegmentNumber, HOMO, Length, COMXCoord, COMYCoord, COMZCoord]
                csvWriter.writerow([segmentNumbers[0], importantEnergies[1], segmentLengths[segmentNumbers[0]], segmentCOM[0], segmentCOM[1], segmentCOM[2]])
                csvFile.close()
            else:
                # print "Segment already present in csv."
                continue


        try:
            if len(failedSingles) > 0:
                print "Calculations completed for single-segments, however there were", len(failedSingles), "errors in calculating the single-segment HOMO levels."
                if len(failedSinglesDict) > 0:
                    print "\n ---======---"
                    print "Failed Singles Dict =", failedSinglesDict
                    print "---======--- \n"
                else:
                    print failedSingles

                jobsToWaitFor = []
                jobsRunning = 0

                for fileName in failedSingles:
                    if fileName not in failedSinglesDict:
                        failedSinglesDict[fileName] = 1
                    elif failedSinglesDict[fileName] == 18:
                        continue
                    else:
                        failedSinglesDict[fileName] += 1
                    failedRerunCounter = failedSinglesDict[fileName]
                    if failedRerunCounter == 3:
                        # Three lots of reruns without any successes, try to turn off SOSCF
                        print str(fileName)+": Three lots of reruns without any success - turning off SOSCF to see if that helps..."
                        turnOffSOSCF(fileName, fullInputPath)
                    if failedRerunCounter == 6:
                        # Still no joy - increase the number of SCF iterations and see if convergence was just slow
                        print str(fileName)+": Six lots of reruns without any success - increasing the number of SCF iterations to 500..."
                        increaseIterations(fileName, fullInputPath)
                    if failedRerunCounter == 9:
                        # Finally, turn down the SCF tolerance
                        print str(fileName)+": Nine lots of reruns without any success - decreasing SCF tolerance (sloppySCF)..."
                        reduceTolerance(fileName, fullInputPath)
                    if failedRerunCounter == 12:
                        print str(fileName)+": Failed to rerun ORCA 12 times, one final thing that can be done is to change the numerical accuracy..."
                        revertORCAFiles(fileName, fullInputPath)
                        increaseGrid(fileName, fullInputPath)
                    if failedRerunCounter == 15:
                        print str(fileName)+": Failed to rerun ORCA 15 times. Will try high numerical accuracy with no SOSCF as a last-ditch effort..."
                        increaseGridNoSOSCF(fileName, fullInputPath)
                    if failedRerunCounter == 18:
                        # SERIOUS PROBLEM
                        print str(fileName)+": Failed to rerun ORCA 18 times, even with all the input file tweaks. Examine the geometry - it is most likely unreasonable."
                        print "Reverting "+str(fileName)+" back to its original state..."
                        revertORCAFiles(fileName, fullInputPath)

                    # NOW RERUN ORCA
                    if failedRerunCounter < 18:
                        if coresToUse == 0:
                            parallel = False
                            # jobNumber, outputFiles, SGEName = configureOrcaSubmission(fileName, fullInputPath+fileName[:-4]+'.inp', orcaPath, outputFiles, engine, queue, parallel)
                            orcaProcess = mp.Process(target=testRun, args=(fileName, fullInputPath+fileName[:-4]+'.inp', orcaPath))
                            orcaProcess.start()
                            orcaProcess.join()
                        else:
                            orcaProcess = mp.Process(target=testRun, args=(fileName, fullInputPath+fileName[:-4]+'.inp', orcaPath))
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
                                            jobsRunning -= 1
                                            jobFinished = True
                                    popList = sorted(popList, reverse=True)
                                    for popIndex in popList:
                                        jobsToWaitFor.pop(popIndex)
                                    if jobFinished == False:
                                        T.sleep(2)


                            # while True:
                            #     parallel = True
                            #     if jobsRunning < coresToUse:
                            #         jobNumber, outputFiles, SGEName = configureOrcaSubmission(fileName, fullInputPath+fileName[:-4]+'.inp', orcaPath, outputFiles, engine, queue, parallel)
                            #         jobsRunning += 1
                            #         jobsToWaitFor.append(jobNumber)
                            #         if engine == 'sge':
                            #             removeSGEs = sgesPresent('./')
                            #             if removeSGEs == True:
                            #                 os.system("rm *.sge.o*")
                            #         break
                            #     else:
                            #         print "\rWaiting for submission slot...",
                            #         T.sleep(2)
                            #         runningJobs, queuedJobs, runningJobNumbers = checkRunningJobs(engine, queue)
                            #         for jobNumber in jobsToWaitFor:
                            #             if jobNumber not in runningJobNumbers:
                            #                 jobsToWaitFor.remove(jobNumber)
                            #                 jobsRunning -= 1
                            # print "\n"

                # Now wait for the jobs to finish
                previousNumberOfJobsToWaitFor = len(jobsToWaitFor)
                while True:
                    popList = []
                    for jobNo in range(len(jobsToWaitFor)):
                        if jobsToWaitFor[jobNo].is_alive() == False:
                            popList.append(jobNo)
                            jobsRunning -= 1
                    popList = sorted(popList, reverse=True)
                    for popIndex in popList:
                        jobsToWaitFor.pop(popIndex)
                    if len(jobsToWaitFor) == 0:
                        break
                    if len(jobsToWaitFor) != previousNumberOfJobsToWaitFor:
                        print "Waiting for", len(jobsToWaitFor), "jobs to finish..."
                        previousNumberOfJobsToWaitFor = len(jobsToWaitFor)
                    T.sleep(2)

                    # runningJobs, queuedJobs, runningJobNumbers = checkRunningJobs(engine, queue)
                    # for jobNumber in jobsToWaitFor:
                    #     if jobNumber not in runningJobNumbers:
                    #         finishedJobs.append(jobNumber)
                    # for jobNumber in finishedJobs:
                    #     jobsToWaitFor.remove(jobNumber)
                    # if len(jobsToWaitFor) == 0:
                    #     break
                    # print "Waiting for", len(jobsToWaitFor), "jobs to finish..."
                    # T.sleep(2)

            # Now check to see if we need to loop through again
            if (len(failedSinglesDict) == 0):
                print "Calculations completed for single-segments with no errors."
                break
            else:
                allowedToBreak = True # Out of the main while loop
                for failCount in failedSinglesDict.values():
                    if failCount != 18: # If ANY of the failed counts in the failedSinglesDict has not reached the maximum (18), then recheck everything.
                        allowedToBreak = False
                        break
                if allowedToBreak == True:
                    break
        except KeyboardInterrupt:
            print "Kill command recieved. Reverting ORCA files..."
            for inputName, failCount in failedSinglesDict.items():
                revertORCAFiles(inputName, fullInputPath)
            print "File reversion complete. Terminating..."
            sys.exit(0)



    # Should now have all of the single-segment files written to a CSV and present in a dictionary
    # Now can start the dual-segment bit


    print "Continuing to double-segment calculations..."


    transferIntegralFileName = fullPath+'TransferIntegrals.csv'

    csvReadOnly = open(transferIntegralFileName, 'a+')
    csvReadOnly.seek(0)
    csvData = csv.reader(csvReadOnly, delimiter = ',')
    for row in csvData:
        HOMOSplitting[(int(row[0]), int(row[1]))] = [float(row[2]), float(row[3]), float(row[4])]
    csvReadOnly.close()

    # failedRerunCounter = 0 # A counter that increments if the number of failed segments did not decrease after one full set of reruns
    # moreIter = False
    # previousLengthOfFailedDoubles = 0


    failedDoublesDict = {} # A dictionary of the names of the failed doubles and the number of failures for each
    outputFiles = None
    koopmanErrors = 0

    while True:
        failedDoubles = []
        for orcaOutput in doubleSegOuts:
            segmentNumbers = obtainSegmentNumbers(orcaOutput[:-4])
            if ((segmentNumbers[0], segmentNumbers[1]) in HOMOSplitting) == False:
                print "Analysing ORCA double-segment file:", orcaOutput
                print "Splitting between segments not present in csv, analysing ORCA file..."
                molecularOrbitals = getOrbitals(fullPath+orcaOutput)
                if molecularOrbitals == False:
                    failedDoubles.append(orcaOutput)
                    continue
                try:
                    # importantEnergies = [segment1ID, segment2ID, HOMO, HOMO-1]
                    importantEnergies = obtainImportantEnergies(segmentNumbers, molecularOrbitals)
                except:
                    failedDoubles.append(orcaOutput)
                    continue
                # Remove the failure state for this output file if it exists
                failedDoublesDict.pop(orcaOutput, None)
                transferIntegral, error = calculateTransferIntegral(segmentNumbers[0], segmentNumbers[1], importantEnergies[2], importantEnergies[3], HOMOEnergies[segmentNumbers[0]], HOMOEnergies[segmentNumbers[1]], koopmansApproximation)
                koopmanErrors += error
                HOMOSplitting[(segmentNumbers[0], segmentNumbers[1])] = [importantEnergies[2], importantEnergies[3], transferIntegral]
                csvFile = open(transferIntegralFileName, 'a')
                csvWriter = csv.writer(csvFile, delimiter = ',')
                csvWriter.writerow([segmentNumbers[0], segmentNumbers[1], importantEnergies[2], importantEnergies[3], transferIntegral])
                csvFile.close()
            else:
                # print "Splitting between these segments already present in csv."
                continue

            # Analyse ORCA
            # Check that the SCF converged
            # Read in the lines and find the MOs
            # Output a CSV with [segment1, segment2, segment1_HOMO, segment1_HOMO-1, segment1_SOLO_HOMO, segment2_HOMO, segment2_HOMO-1, segment2_SOLO_HOMO]
            # Or maybe just [segment1, segment2, transferIntegral]
            # This can then be read in by chargeTransport.py

            # If len(importantEnergies) == 2: Put it in a dict - HOMO levels only
            # if len == 4: Keep in list, iterate over this list to generate the TIs and put them in the CSV


        # print "TOTAL NUMBER OF KOOPMAN'S ERRORS =", koopmanErrors

        try:
            if len(failedDoubles) > 0:
                print "\nCalculations completed for double-segments, however there were", len(failedDoubles), "errors in calculating the double-segment HOMO splitting."
                if len(failedDoublesDict) > 0:
                    print "\n ---======---"
                    print "Failed Doubles Dict =", failedDoublesDict
                    print "---======--- \n"
                else:
                    print failedDoubles


                jobsToWaitFor = []
                jobsRunning = 0


                # # TEMPORARY
                # for fileName in failedDoubles:
                #     revertORCAFiles(fileName, fullInputPath)
                # print "Files successfully reverted."
                # sys.exit(0)




                for fileName in failedDoubles:
                    if fileName not in failedDoublesDict:
                        failedDoublesDict[fileName] = 1
                    elif failedDoublesDict[fileName] == 18:
                        continue
                    else:
                        failedDoublesDict[fileName] += 1
                    failedRerunCounter = failedDoublesDict[fileName]
                    if failedRerunCounter == 3:
                        # Three lots of reruns without any successes, try to turn off SOSCF
                        print str(fileName)+": Three lots of reruns without any success - turning off SOSCF to see if that helps..."
                        turnOffSOSCF(fileName, fullInputPath)
                    if failedRerunCounter == 6:
                        # Still no joy - increase the number of SCF iterations and see if convergence was just slow
                        print str(fileName)+": Six lots of reruns without any success - increasing the number of SCF iterations to 500..."
                        increaseIterations(fileName, fullInputPath)
                    if failedRerunCounter == 9:
                        # Finally, turn down the SCF tolerance
                        print str(fileName)+": Nine lots of reruns without any success - decreasing SCF tolerance (sloppySCF)..."
                        reduceTolerance(fileName, fullInputPath)
                    if failedRerunCounter == 12:
                        print str(fileName)+": Failed to rerun ORCA 12 times, one final thing that can be done is to change the numerical accuracy..."
                        revertORCAFiles(fileName, fullInputPath)
                        increaseGrid(fileName, fullInputPath)
                    if failedRerunCounter == 15:
                        print str(fileName)+": Failed to rerun ORCA 15 times. Will try high numerical accuracy with no SOSCF as a last-ditch effort..."
                        increaseGridNoSOSCF(fileName, fullInputPath)
                    if failedRerunCounter == 18:
                        # SERIOUS PROBLEM
                        print str(fileName)+": Failed to rerun ORCA 18 times, even with all the input file tweaks. Examine the geometry - it is most likely unreasonable."
                        print "Reverting "+str(fileName)+" back to its original state..."
                        revertORCAFiles(fileName, fullInputPath)

                    # NOW RERUN ORCA
                    if failedRerunCounter < 18:
                        if coresToUse == 0:
                            parallel = False
                            # jobNumber, outputFiles, SGEName = configureOrcaSubmission(fileName, fullInputPath+fileName[:-4]+'.inp', orcaPath, outputFiles, engine, queue, parallel)
                            orcaProcess = mp.Process(target=testRun, args=(fileName, fullInputPath+fileName[:-4]+'.inp', orcaPath))
                            orcaProcess.start()
                            orcaProcess.join()
                        else:
                            orcaProcess = mp.Process(target=testRun, args=(fileName, fullInputPath+fileName[:-4]+'.inp', orcaPath))
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
                                            jobsRunning -= 1
                                            jobFinished = True
                                    popList = sorted(popList, reverse=True)
                                    for popIndex in popList:
                                        jobsToWaitFor.pop(popIndex)
                                    if jobFinished == False:
                                        T.sleep(2)

                            # printWaiting = False
                            # while True:
                            #     parallel = True
                            #     if jobsRunning < coresToUse:
                            #         jobNumber, outputFiles, SGEName = configureOrcaSubmission(fileName, fullInputPath+fileName[:-4]+'.inp', orcaPath, outputFiles, engine, queue, parallel)
                            #         jobsRunning += 1
                            #         # print "New Job submitted: JobsRunning ==", jobsRunning, "coresToUse ==", coresToUse
                            #         jobsToWaitFor.append(jobNumber)
                            #         if engine == 'sge':
                            #             removeSGEs = sgesPresent('./')
                            #             if removeSGEs == True:
                            #                 os.system("rm *.sge.o*")
                            #         printWaiting = True
                            #         break
                            #     else:
                            #         if (printWaiting == True):
                            #             printWaiting = False
                            #             print "Waiting for submission slot..."# Jobs to Wait for =", jobsToWaitFor, "RunningJobNumbers =", runningJobNumbers
                            #         T.sleep(2)
                            #         runningJobs, queuedJobs, runningJobNumbers = checkRunningJobs(engine, queue)
                            #         for jobNumber in jobsToWaitFor:
                            #             if jobNumber not in runningJobNumbers:
                            #                 jobsToWaitFor.remove(jobNumber)
                            #                 jobsRunning -= 1
                            # # print "\n"

                # Now wait for the jobs to finish
                previousNumberOfJobsToWaitFor = len(jobsToWaitFor)
                while True:
                    popList = []
                    for jobNo in range(len(jobsToWaitFor)):
                        if jobsToWaitFor[jobNo].is_alive() == False:
                            popList.append(jobNo)
                            jobsRunning -= 1
                    popList = sorted(popList, reverse=True)
                    for popIndex in popList:
                        jobsToWaitFor.pop(popIndex)
                    if len(jobsToWaitFor) == 0:
                        break
                    if len(jobsToWaitFor) != previousNumberOfJobsToWaitFor:
                        print "Waiting for", len(jobsToWaitFor), "jobs to finish..."
                        previousNumberOfJobsToWaitFor = len(jobsToWaitFor)
                    T.sleep(2)


                # previousLenJobs = 0
                # while True:
                #     finishedJobs = []
                #     runningJobs, queuedJobs, runningJobNumbers = checkRunningJobs(engine, queue)
                #     for jobNumber in jobsToWaitFor:
                #         if jobNumber not in runningJobNumbers:
                #             finishedJobs.append(jobNumber)
                #     for jobNumber in finishedJobs:
                #         jobsToWaitFor.remove(jobNumber)
                #     if len(jobsToWaitFor) == 0:
                #         break
                #     if (len(jobsToWaitFor) != previousLenJobs):
                #         previousLenJobs = len(jobsToWaitFor)
                #         print "Waiting for", len(jobsToWaitFor), "jobs to finish..."
                #     T.sleep(2)


            if (len(failedDoublesDict) == 0):
                print "Calculations completed for double-segments with no errors."
                break
            else:
                allowedToBreak = True # Out of main while loop
                for failCount in failedDoublesDict.values():
                    if failCount != 18:
                        allowedToBreak = False
                        break # Out of for loop
                if allowedToBreak == True:
                    break
        except KeyboardInterrupt:
            print "Kill command recieved. Reverting ORCA files..."
            for inputName, failCount in failedDoublesDict.items():
                if (failCount != 18) and (failCount > 3):
                    revertORCAFiles(inputName, fullInputPath)
            print "File reversion complete. Terminating..."
            sys.exit(0)



    if len(failedSingles) != 0:
        print "All", len(failedSingles), "final single segment failures for", morphology
        print sorted(failedSinglesDict.keys())
    if len(failedDoubles) != 0:
        print "All", len(failedDoubles), "final double segment failures for", morphology
        print sorted(failedDoublesDict.keys())
