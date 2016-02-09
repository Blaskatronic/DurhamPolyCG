import os
import numpy as np
import pylab as P
import time as T
import csv
import subprocess
import argparse
import copy

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


def turnOffSOSCF(failedFiles, fullInputPath):
    for outputFile in failedFiles:
        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
        originalLines = fileName.readlines()
        fileName.close()
        
        originalLines[3] = '!ZINDO/S NoSOSCF\n'

        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
        fileName.writelines(originalLines)
        fileName.close()


def reduceTolerance(failedFiles, fullInputPath):
    for outputFile in failedFiles:
        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
        originalLines = fileName.readlines()
        fileName.close()

        originalLines[3] = '!ZINDO/S NoSOSCF SloppySCF\n'

        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
        fileName.writelines(originalLines)
        fileName.close()


def increaseIterations(failedFiles, fullInputPath):
    for outputFile in failedFiles:
        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
        originalLines = fileName.readlines()
        fileName.close()

        originalLines.append('\n%scf MaxIter 500 end')

        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
        fileName.writelines(originalLines)
        fileName.close()


def increaseGrid(failedFiles, fullInputPath):
    for outputFile in failedFiles:
        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
        originalLines = fileName.readlines()
        fileName.close()

        originalLines[3] = '!ZINDO/S SlowConv Grid7 NoFinalGrid\n'

        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
        fileName.writelines(originalLines)
        fileName.close()


def increaseGridNoSOSCF(failedFiles, fullInputPath):
    for outputFile in failedFiles:
        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
        originalLines = fileName.readlines()
        fileName.close()
        
        originalLines[3] = '!ZINDO/S SlowConv Grid7 NoFinalGrid NoSOSCF SloppySCF\n'
        originalLines.append('\n%scf MaxIter 500 end')

        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
        fileName.writelines(originalLines)
        fileName.close()


def revertORCAFiles(files, fullInputPath, moreIter):
    for outputFile in files:
        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'r')
        originalLines = fileName.readlines()
        fileName.close()

        originalLines[3] = '! ZINDO/S\n'
        if moreIter == True:
            originalLines.pop(-1)

        fileName = open(fullInputPath+outputFile[:-3]+'inp', 'w+')
        fileName.writelines(originalLines)
        fileName.close()


def rerunOrca(failedSingles, orcaPath, inputPath, outputPath):
    for outputFile in failedSingles:
        runOrca = orcaPath+'/orca '+inputPath+outputFile[:-3]+'inp > '+outputPath+outputFile
        print runOrca
        subprocess.call(runOrca, shell=True)
    


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


    failedRerunCounter = 0
    moreIter = False
    previousLengthOfFailedSingles = 0
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
                importantEnergies = obtainImportantEnergies(segmentNumbers, molecularOrbitals)
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



        if len(failedSingles) > 0:
            print "Calculations completed for single-segments, however there were", len(failedSingles), "errors in calculating the single-segment HOMO levels."
            print "The files that failed were:"
            print failedSingles
            print "Rerunning failed Singles:"
            if len(failedSingles) == previousLengthOfFailedSingles:
                failedRerunCounter += 1
                print "Failed Rerun counter =", failedRerunCounter
            else:
                revertORCAFiles(failedSingles, fullInputPath, moreIter)
                moreIter = False
                failedRerunCounter = 0
            if failedRerunCounter == 3:
                # Three lots of reruns without any successes, try to turn off SOSCF
                print "Three lots of reruns without any success - turning off SOSCF to see if that helps..."
                turnOffSOSCF(failedSingles, fullInputPath)
            if failedRerunCounter == 6:
                # Still no joy - increase the number of SCF iterations and see if convergence was just slow
                print "Six lots of reruns without any success - increasing the number of SCF iterations to 500..."
                moreIter = True
                increaseIterations(failedSingles, fullInputPath)
            if failedRerunCounter == 9:
                # Finally, turn down the SCF tolerance
                print "Nine lots of reruns without any success - decreasing SCF tolerance (sloppySCF)..."
                reduceTolerance(failedSingles, fullInputPath)
            if failedRerunCounter == 12:
                print "Failed to rerun ORCA 12 times, one final thing that can be done is to change the numerical accuracy..."
                revertORCAFiles(failedSingles, fullInputPath)
                increaseGrid(failedSingles, fullInputPath)
                moreIter = False
                pass
            if failedRerunCounter == 15:
                print "Failed to rerun ORCA 15 times. Will try high numerical accuracy with no SOSCF as a last-ditch effort..."
                increaseGridNoSOSCF(failedSingles, fullInputPath)
                moreIter = True
            if failedRerunCounter == 18:
                # SERIOUS PROBLEM
                print "Failed to rerun ORCA 18 times, even with all the input file tweaks. Examine the geometry - it is most likely unreasonable."
                print "Reverting the input files back to their original state..."
                revertORCAFiles(failedSingles, fullInputPath, moreIter)
                print "Continuing to double-segment calculations..."
                break
            rerunOrca(failedSingles, orcaPath, fullInputPath, fullPath)

        else:
            print "Calculations completed for single-segments with no errors."
            break


    # Should now have all of the single-segment files written to a CSV and present in a dictionary




    # Now can start the dual-segment bit

    transferIntegralFileName = fullPath+'TransferIntegrals.csv'

    csvReadOnly = open(transferIntegralFileName, 'a+')
    csvReadOnly.seek(0)
    csvData = csv.reader(csvReadOnly, delimiter = ',')
    for row in csvData:
        HOMOSplitting[(int(row[0]), int(row[1]))] = [float(row[2]), float(row[3]), float(row[4])]
    csvReadOnly.close()

    failedRerunCounter = 0 # A counter that increments if the number of failed segments did not decrease after one full set of reruns
    moreIter = False
    previousLengthOfFailedDoubles = 0
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
                    importantEnergies = obtainImportantEnergies(segmentNumbers, molecularOrbitals)
                except:
                    failedDoubles.append(orcaOutput)
                    continue
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

        if len(failedDoubles) > 0:
            print "Calculations completed for double-segments, however there were", len(failedDoubles), "errors in calculating the double-segment HOMO splitting."
            print "The files that failed were:"
            print failedDoubles
            print "Rerunning failed Doubles:"
            if len(failedDoubles) == previousLengthOfFailedDoubles:
                failedRerunCounter += 1
                print "Failed Rerun counter =", failedRerunCounter
            else:
                revertORCAFiles(failedDoubles, fullInputPath, moreIter)
                moreIter = False
                failedRerunCounter = 0
            if failedRerunCounter == 3:
                # Three lots of reruns without any successes, try to turn off SOSCF
                print "Three lots of reruns without any success - turning off SOSCF to see if that helps..."
                turnOffSOSCF(failedDoubles, fullInputPath)
            if failedRerunCounter == 6:
                # Still no joy - increase the number of SCF iterations and see if convergence was just slow
                print "Six lots of reruns without any success - increasing the number of SCF iterations to 500..."
                moreIter = True
                increaseIterations(failedDoubles, fullInputPath)
            if failedRerunCounter == 9:
                # Finally, turn down the SCF tolerance
                print "Nine lots of reruns without any success - decreasing SCF tolerance (sloppySCF)..."
                reduceTolerance(failedDoubles, fullInputPath)
            if failedRerunCounter == 12:
                print "Failed to rerun ORCA 12 times, one final thing that can be done is to change the numerical accuracy..."
                revertORCAFiles(failedDoubles, fullInputPath, moreIter)
                moreIter = False
                increaseGrid(failedDoubles, fullInputPath)
            if failedRerunCounter == 15:
                print "Failed to rerun ORCA 15 times. Will try high numerical accuracy with no SOSCF as a last-ditch effort..."
                moreIter = True
                increaseGridNoSOSCF(failedDoubles, fullInputPath)
            if failedRerunCounter == 18:
                # SERIOUS PROBLEM
                print "Failed to rerun ORCA 18 times, even with all the input file tweaks. Examine the geometry - it is most likely unreasonable."
                print "Reverting the input files back to their original state..."
                revertORCAFiles(failedDoubles, fullInputPath, moreIter)
                print "Done. Exiting program."
                break
            rerunOrca(failedDoubles, orcaPath, fullInputPath, fullPath)


            previousLengthOfFailedDoubles = len(failedDoubles)

        else:
            print "All calculations completed with no errors."
            break

    if len(failedSingles) != 0:
        print "Failed Singles =", failedSingles
    if len(failedDoubles) != 0:
        print "Failed Doubles =", failedDoubles
