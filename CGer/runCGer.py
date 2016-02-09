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
import sys


def getFilesList(direc, mode='both'):
    fileList = os.listdir(direc)
#    files = []
#    for element in fileList:
#        files.append(element[7:])
    datFiles = []
    lammpsFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.dat'):
            datFiles.append(str(fileName))

    if mode == 'both':
        popList = []
        for datNumber in range(len(datFiles)):
            for fileName in fileList:
                addedLammps = 0
                if (datFiles[datNumber] in fileName) and ('.lammpstrj' in fileName):
                    lammpsFiles.append(fileName)
                    addedLammps = 1
                    break
            if addedLammps == 0:
                popList.append(datNumber)
        popList = sorted(popList, reverse=True)
        for popIndex in popList:
            datFiles.pop(popIndex)

        validFiles = sorted(zip(datFiles, lammpsFiles))

    else:
        validFiles = sorted(datFiles)
    
    return validFiles


def modifiedSchulzFlory(lengthBase, lengthFactor):
    multiplicityFactor = (R.random()**lengthFactor)*lengthFactor
    chainMultiplied = lengthBase*multiplicityFactor
    monomersInChain = 2*(np.round(chainMultiplied/2.))
    if monomersInChain == 0:
        monomersInChain = 2
    return monomersInChain


def calcSFDistribution(tunableParameter):
    sfDist = []
    sfCumDist = []
    runningTotal = 0
    monomerLengths = list(np.arange(0,401,1)) # Need to do each individual step to get the right cumulative distribution
    for length in monomerLengths:
        probability = tunableParameter**2 * length * ((1 - tunableParameter) ** (length - 1))
        sfDist.append(probability)
        runningTotal += probability
        sfCumDist.append(runningTotal)

    sfCumDist.pop(0) # The one representing a 0 length chain
    monomerLengths.pop(0)
    sfCumDist.pop(0) # The one representing a 1 length chain
    monomerLengths.pop(0)
    popList = list(np.arange(1, 399, 2))
    popList.sort(reverse = True)
    for i in popList:
        sfCumDist.pop(i)
        monomerLengths.pop(i)
    return sfCumDist


def unmodifiedSchulzFlory(SFDistribution, lengthBase):
    chainLengths = list(np.arange(2,lengthBase,2))
    randomSelection = R.random()
    for elementNo in range(len(SFDistribution)):
        if randomSelection < SFDistribution[elementNo]:
            monomerIndex = elementNo
            break
        else:
            monomerIndex = elementNo
    return chainLengths[monomerIndex]


def calculateP3HTChainParameters(chainLengths):
    massOfOneMonomer = 166.28708 #au for P3HT
    totalChainMass = 0
    totalChainMass_2 = 0

    for chainLength in chainLengths:
        chainMass = massOfOneMonomer*chainLength
        totalChainMass += chainMass
        totalChainMass_2 += chainMass**2

    Mn = totalChainMass/float(numberOfChains)
    Mw = totalChainMass_2/totalChainMass
    PDI = Mw/Mn

    return Mn, Mw, PDI


def determineMixture(direc):
    allSampleFiles = getFilesList(direc, mode='datOnly')
    sampleFiles = copy.deepcopy(allSampleFiles)
    while True:
        print "\n---=== SAMPLES ON THE SAMPLE SHELF THAT CAN BE MIXED ===---"
        if len(sampleFiles) == 0:
            print "ERROR: No valid samples found in", str(direc)+". Please check that the dat files for the samples to be mixed are present in this folder before running!"
            return None
        for sampleNo in range(len(sampleFiles)):
            print str(sampleNo)+"):", sampleFiles[sampleNo]
        print str(sampleNo+1)+"): Exit CGMixer"
        mixOneNo = raw_input("Please select the first sample to be mixed (integer, default = 0): ")
        if len(mixOneNo) == 0:
            mixOneNo = 0
        else:
            try:
                mixOneNo = int(mixOneNo)
            except:
                print "Please enter an integer between 0 and", len(sampleFiles)
                continue
        if (mixOneNo < 0) or (mixOneNo > len(sampleFiles)):
            print "Please enter an integer between 0 and", len(sampleFiles)
            continue
        elif mixOneNo == len(sampleFiles):
            print "Exiting..."
            return None
        else:
            break
    while True:
        prop1 = raw_input("Please indicate the proportion of this sample that you would like to use in the final mixture (float, default = 1.0): ")
        print "---=======================================================---"
        if len(prop1) == 0:
            prop1 = 1.0
        else:
            try:
                prop1 = float(prop1)
            except:
                print "Please enter a float between 0.0 and 1.0"
                continue
        if (prop1 < 0.0) or (prop1 > 1.0):
            print "Please enter a float between 0.0 and 1.0"
            continue
        elif (prop1 == 1.0):
            return [allSampleFiles[mixOneNo], prop1, 0, 0]
        else:
            break
    sampleFiles.pop(mixOneNo)
    while True:
        print "\n---=== SAMPLES ON THE SAMPLE SHELF THAT CAN BE MIXED ===---"
        if len(sampleFiles) == 0:
            print "ERROR: No additional samples found in", str(direc)+"! Please check that the desired dat files for all samples are present in this folder before running!"
            return [allSampleFiles[mixOneNo], prop1, 0, 0]
        for sampleNo in range(len(sampleFiles)):
            print str(sampleNo)+"):", sampleFiles[sampleNo]
        print str(sampleNo+1)+"): Do not mix a second sample"
        mixTwoNo = raw_input("Please select the second sample to be mixed (integer, default = 0): ")
        if len(mixTwoNo) == 0:
            mixTwoNo = 0
        else:
            try:
                mixTwoNo = int(mixTwoNo)
            except:
                print "Please enter an integer between 0 and", len(sampleFiles)
                continue
        if (mixTwoNo < 0) or (mixTwoNo > len(sampleFiles)):
            print "Please enter an integer between 0 and", len(sampleFiles)
            continue
        elif mixTwoNo == len(sampleFiles):
            return [allSampleFiles[mixOneNo], prop1, 0, 0]
        else:
            break
    # print "This two-sample mixture has", str(prop1*100)+"% of sample 1 and so therefore will have", str((1-prop1)*100)+"% of sample 2."
    print "---=======================================================---"
    return [allSampleFiles[mixOneNo], prop1, sampleFiles[mixTwoNo], 1-prop1]


def generateMixture(chainLengthsDict, mass1, prop1, mass2, prop2):
    newMixture = []
    if mass1 != None:
        sample1 = chainLengthsDict[mass1]
        numberOfM1Chains = np.round(prop1*len(sample1))
        for i in range(int(numberOfM1Chains)):
            chainLength = R.choice(sample1)
            sample1.remove(chainLength)
            newMixture.append(chainLength)
    if mass2 != None:
        sample2 = chainLengthsDict[mass2]
        numberOfM2Chains = np.round(prop2*len(sample2))
        for i in range(int(numberOfM2Chains)):
            chainLength = R.choice(sample2)
            sample2.remove(chainLength)
            newMixture.append(chainLength)

    return newMixture



def loadDat(filename):
    print "Reading in .dat file:", filename
    bonds = []
    angles = []

    originalFile = open('./'+str(filename), 'r')
    datData = originalFile.readlines()
    originalFile.close()

    header = datData[1:7]

    atomNo, bondNo, angleNo, dihedralNo, improperNo = treatHeader(header)

    for lineNo in range(len(datData)):
        if "Masses" in datData[lineNo]:
            massStart = lineNo+2
        if "Atoms" in datData[lineNo]:
            massEnd = lineNo-1
            atomStart = lineNo+2
        if "Bonds" in datData[lineNo]:
            bondStart = lineNo+2
        if "Angles" in datData[lineNo]:
            angleStart = lineNo+2
        if "Dihedrals" in datData[lineNo]:
            dihedralStart = lineNo+2
        if "Impropers" in datData[lineNo]:
            improperStart = lineNo+2

    atomEnd = atomStart+atomNo
    bondEnd = bondStart+bondNo
    angleEnd = angleStart+angleNo
    dihedralEnd = dihedralStart+dihedralNo
    improperEnd = improperStart+improperNo

    masses = datData[massStart:massEnd]
    masses = treatDat(masses, 'masses')
    atoms = datData[atomStart:atomEnd]
    atoms = treatDat(atoms, 'atoms')
    bonds = datData[bondStart:bondEnd]
    bonds = treatDat(bonds, 'bonds')
    angles = datData[angleStart:angleEnd]
    angles = treatDat(angles, 'angles')
    dihedrals = datData[dihedralStart:dihedralEnd]
    dihedrals = treatDat(dihedrals, 'dihedrals')
    impropers = datData[improperStart:improperEnd]
    impropers = treatDat(impropers, 'impropers')

    return masses, atoms, bonds, angles, dihedrals, impropers

def treatDat(data, dataType):
    newData = []
#    print data[0]
    if dataType == 'masses':
        for element in data:
            newData.append([])
            tempLine = element[:-1].split(' ')
            for splitElement in tempLine:
                if len(splitElement) != 0:
                    newData[-1].append(float(splitElement))
        for massNo in range(len(newData)):
            newData[massNo][0] = int(newData[massNo][0])
    elif dataType == 'atoms':
        for element in data:
            newData.append([])
            tempLine = element[:-1].split(' ')
            for splitElement in tempLine:
                if len(splitElement) != 0:
                    newData[-1].append(splitElement)

        for atomNo in range(len(newData)):
            for valueNo in range(len(newData[atomNo])):
                if (valueNo == 3) or (valueNo == 4) or (valueNo == 5):
                    if "e" in newData[atomNo][valueNo]:
                        newData[atomNo][valueNo] = float(0)
                    else:
                        newData[atomNo][valueNo] = float(newData[atomNo][valueNo])
                else:
                    newData[atomNo][valueNo] = int(newData[atomNo][valueNo])
    else:
        for element in data:
            newData.append([])
            tempLine = element[:-1].split(' ')
            for splitElement in tempLine:
                if len(splitElement) != 0:
                    newData[-1].append(int(splitElement))

    return newData


def treatHeader(header):
    numbersList = []
    for line in header:
        newLine = []
        tempLine = line.split(' ')
        for element in tempLine:
            try:
                numbersList.append(int(element))
            except:
                continue
    return numbersList[0], numbersList[1], numbersList[2], numbersList[3], numbersList[4]


def getChainLengths(atoms):
    moleculeLengths = []
    numberOfCGAtomsInMolecule = {}
    for atom in atoms:
        chainNumber = atom[1]
        if chainNumber not in numberOfCGAtomsInMolecule:
            numberOfCGAtomsInMolecule[chainNumber] = 0
        numberOfCGAtomsInMolecule[chainNumber] += 1
    # MoleculeLengths now have the number of CG atoms within the chain. Need to get this down to number of monomers.
    for moleculeNo in numberOfCGAtomsInMolecule:
        moleculeLength = numberOfCGAtomsInMolecule[moleculeNo]/3
        moleculeLengths.append(moleculeLength)
    return moleculeLengths
    
        

def treatDatName(datName):
    underscoreLocs = findIndex(datName, '_')
    Mw = datName[underscoreLocs[0]+1:underscoreLocs[1]]
    return Mw


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
    exitFlag = 0
    while exitFlag == 0:
        while True:
            print "\n---=== MODES THAT CAN BE EXECUTED BY CGER.PY ===---"
            print "0) Create initial sample of chains (using modified Schulz-Flory distribution) to be submitted to LAMMPS for 'fastshrinking'"
            print "1) Load 'fastshrunked' morphology to be resubmitted into LAMMPS for 'heating' process"
            print "2) Load 'heated' morphology to be resubmitted into LAMMPS for 'cooling' process"
            print "3) Load 'cooled' morphology to be resubmitted into LAMMPS for 'equilibration' process"
            print "4) Create Single Chain of Specified Length (for lookup table)"
            print "5) Create an mixed sample of chains from two previously generated samples to be submitted to LAMMPS for 'fastshrinking'"
            print "6) Exit program"
            mode = raw_input("Please select an execution mode (integer, default = 0): ")
            print "---===============================================---"
            if len(mode) == 0:
                mode = 0
                break
            else:
                try:
                    mode = int(mode)
                    break
                except:
                    print "Please enter an integer between 0 and 6."
                    continue
            if (mode < 0) or (mode > 6):
                print "Please enter an integer between 0 and 6."
                continue
            elif mode == 6:
                print "Exiting program..."
                exitFlag = 1
                break
        if (mode == 0): # Create normal morphology
            while True:
                # Set default parameters:
                MSF = True
                lengthBase = 280
                lengthFactor = 1.6
                # lengthBase = 200
                # lengthFactor = 1.9
                schulzFloryIndex = 0.008
                numberOfChains = 40
                chainSpecies = 'con'
                chainRotation = True
                makeLammps = True
                plotGraph = True
                chainLengths = []
                if MSF == False:
                    SFDist = calcSFDistribution(schulzFloryIndex)
                for i in range(numberOfChains):
                    if MSF == True:
                        monomersInChain = modifiedSchulzFlory(lengthBase, lengthFactor)
                    else:
                        monomersInChain = unmodifiedSchulzFlory(SFDist, lengthBase)
                    chainLengths.append(int(monomersInChain))
                Mn, Mw, PDI = calculateP3HTChainParameters(chainLengths)
                print "\n---=== Current Chain Parameters: ===---"
                print "Use Modified Schulz-Flory Distributon (MSF) =", MSF
                print "Length Base (maximum chain length when MSF = False, approximate chain length when MSF = True) =", lengthBase
                print "Length Factor (Unused when MSF = False, helps calculate chain deviation when MSF = True) =", lengthFactor
                print "Schulz-Flory Index (Used to calculate SF distribution when MSF = False) =", schulzFloryIndex
                print "Resultant Mn =", Mn
                print "Resultant Mw =", Mw
                print "Resultant PDI =", PDI
                print "Number of Chains in Sample =", numberOfChains
                print "Chain Species =", chainSpecies
                print "Chain rotation =", chainRotation
                print "Make LAMMPS input .dat file =", makeLammps
                print "Plot a 3D graph of the atoms in the resultant .dat file =", plotGraph
                print "---==================================---"

                while True:
                    print "\n0) Accept current parameters and continue (default)"
                    print "1) Recalculate this chain length distribution"
                    paramOption = raw_input("Please choose whether to modify parameters (integer): ")
                    if (len(paramOption) == 0):
                        paramOption = 0
                        break
                    else:
                        try:
                            paramOption = int(paramOption)
                        except:
                            paramOption = 0
                            break
                    if (paramOption < 0) or (paramOption > 1):
                        paramOption = 0
                        break
                    else:
                        break
                if paramOption == 0:
                    print "Accepting these parameters..."
                    break
                elif paramOption == 1:
                    continue
            # Now we have a bunch of parameters and a "Custom Length" style array of chain lengths.
            t0 = T.time()
            print 'python CGer.py -l "'+str(chainLengths)+'" -s '+str(chainSpecies)+' -r '+str(int(chainRotation))+' -x '+str(mode)+' -d '+str(int(makeLammps))+' -g '+str(int(plotGraph))
            subprocess.call('python CGer.py -l "'+str(chainLengths)+'" -s '+str(chainSpecies)+' -r '+str(int(chainRotation))+' -x '+str(mode)+' -d '+str(int(makeLammps))+' -g '+str(int(plotGraph)), shell=True)
            # Return to main menu
            continue
        elif (mode == 1): # Load fastshrunk
            direc = './shrunkVols/'
            validFiles = getFilesList(direc)
        elif (mode == 2): # Load heated
            direc = './preCooling/'
            validFiles = getFilesList(direc)
        elif (mode == 3): # Load cooled
            direc = './preEq/'
            validFiles = getFilesList(direc)
        elif (mode == 4): # Create Single
            lengthBase = 280
            chainSpecies = 'str'
            chainRotation = False
            makeLammps = True
            plotGraph = True
            print "\n---=== Current Chain Parameters: ===---"
            print "Single Chain"
            print "Chain Species =", chainSpecies
            print "Chain Length =", lengthBase
            print "Chain Rotation =", chainRotation
            print "Make LAMMPS Input .dat file =", makeLammps
            print "Plot a 3D graph of the atoms in the resultant .dat file =", plotGraph
            print "---==================================---"
            while True:
                print "\n0) Accept current parameters and continue (default)"
                print "1) Modify parameters (not working currently)"
                paramOption = raw_input("Please choose whether to modify parameters: ")
                if (len(paramOption) == 0):
                    paramOption = 0
                else:
                    try:
                        paramOption = int(paramOption)
                    except:
                        print "Accepting these parameters..."
                        paramOption = 0
                if (paramOption < 0) or (paramOption > 1):
                    print "Accepting these parameters..."
                    paramOption = 0
                break
            if paramOption == 1:
                raise SystemError('FEATURE NOT ADDED YET')
            else:
                print 'python CGer.py -l "'+str([lengthBase])+'" -s '+str(chainSpecies)+' -r '+str(int(chainRotation))+' -x '+str(mode)+' -d '+str(int(makeLammps))+' -g '+str(int(plotGraph))
                subprocess.call('python CGer.py -l "'+str([lengthBase])+'" -s '+str(chainSpecies)+' -r '+str(int(chainRotation))+' -x '+str(mode)+' -d '+str(int(makeLammps))+' -g '+str(int(plotGraph)), shell=True)
            # Return to main menu
            continue
        elif (mode == 5): # Create Mixed Sample
            direc = './sampleShelf/'
            getMixtureParameters = True
            while True:
                if getMixtureParameters == True:
                    mixtureParameters = determineMixture(direc)
                print "\n"
                if mixtureParameters == None:
                    continue
                else:
                    if (mixtureParameters[1] != 0.0) and (mixtureParameters[3] != 0.0):
                        print "Mixing a sample with", str(mixtureParameters[1]*100)+"% of", str(mixtureParameters[0])+", and", str(mixtureParameters[3]*100)+"% of", str(mixtureParameters[2])+"..."
                        mixedDats = [mixtureParameters[0], mixtureParameters[2]]
                    elif (mixtureParameters[1] == 0.0):
                        print "Mixing a sample with", str(mixtureParameters[3]*100)+"% of", str(mixtureParameters[2])+"..."
                        mixedDats = [mixtureParameters[2]]
                    elif (mixtureParameters[3] == 0.0):
                        print "Mixing a sample with", str(mixtureParameters[1]*100)+"% of", str(mixtureParameters[0])+"..."
                        mixedDats = [mixtureParameters[0]]
                numberOfChains = 40
                chainSpecies = 'con'
                chainRotation = True
                makeLammps = True
                plotGraph = True
                chainLengthsDict = {}
                Mws = [None, None]
                for datNo in range(len(mixedDats)):
                    Mw = treatDatName(mixedDats[datNo])
                    Mws[datNo] = int(Mw)
                    masses, atoms, bonds, angles, dihedrals, impropers = loadDat(direc+mixedDats[datNo])
                    chainLengths = getChainLengths(atoms)
                    chainLengthsDict[int(Mw)] = chainLengths
                customMixtureOfLengths = generateMixture(chainLengthsDict, Mws[0], mixtureParameters[1], Mws[1], mixtureParameters[3])
                Mn, Mw, PDI = calculateP3HTChainParameters(customMixtureOfLengths)
                print "Custom lenghts =", customMixtureOfLengths

                print "\n---=== Current Mixture Parameters: ===---"
                print "Sample 1 =", mixtureParameters[0]
                print "Proportion of Sample 1 =", str(mixtureParameters[1]*100)+"%"
                print "Sample 2 =", mixtureParameters[2]
                print "Proportion of Sample 2 =", str(mixtureParameters[3]*100)+"%"
                print "Resultant Mn =", Mn
                print "Resultant Mw =", Mw
                print "Resultant PDI =", PDI
                print "Number of Chains in Sample =", numberOfChains
                print "Chain Species =", chainSpecies
                print "Chain rotation =", chainRotation
                print "Make LAMMPS input .dat file =", makeLammps
                print "Plot a 3D graph of the atoms in the resultant .dat file =", plotGraph
                print "---==================================---"
                while True:
                    print "\n0) Accept current parameters and continue (default)"
                    print "1) Remix these proportions for more desireable parameters"
                    print "2) Remix a different combination/differet samples"
                    paramOption = raw_input("Please choose whether to modify parameters: ")
                    if (len(paramOption) == 0):
                        paramOption = 0
                        break
                    else:
                        try:
                            paramOption = int(paramOption)
                            break
                        except:
                            paramOption = 0
                            break
                    if (paramOption < 0) or (paramOption > 2):
                        paramOption = 0
                        break
                if paramOption == 0:
                    print "Accepting these parameters..."
                    break
                elif paramOption == 1:
                    print "Rerolling mixture..."
                    getMixtureParameters = False
                    continue
                elif paramOption == 2:
                    print "Creatign new mixture..."
                    getMixtureParamters = True
                    continue
 
            print 'python CGer.py -l "'+str(customMixtureOfLengths)+'" -s '+str(chainSpecies)+' -r '+str(int(chainRotation))+' -x '+str(mode)+' -d '+str(int(makeLammps))+' -g '+str(int(plotGraph))
            subprocess.call('python CGer.py -l "'+str(customMixtureOfLengths)+'" -s '+str(chainSpecies)+' -r '+str(int(chainRotation))+' -x '+str(mode)+' -d '+str(int(makeLammps))+' -g '+str(int(plotGraph)), shell=True)
            print "Calculations complete."
            continue
        else: # Exit
            sys.exit()


        # At this point we are loading a morphology, so can use the loading code

        while True:
            print "\n---=== FILES THAT CAN BE RUN BY CGER.PY ===---"
            if len(validFiles) == 0:
                print "ERROR: No valid file pairs found in", str(direc)+". Please check that both the .lammpstrj and the corresponding .dat are present in this folder before running!"
                exitFlag = 1
                break
            for elementNo in range(len(validFiles)):
                print str(elementNo)+"):", validFiles[elementNo]
            print str(elementNo+1)+"): Exit runCGer.py"
            # print "Valid files =", zip(datFiles, lammpstrjFiles)
            runThisFile = raw_input("Please pick a file to run (integer, default = 0): ")
            if len(runThisFile) == 0:
                runThisFile = 0
            else:
                try:
                    runThisFile = int(runThisFile)
                except:
                    print "Please enter an integer between 0 and", len(validFiles)
                    continue
            if (runThisFile < 0) or (runThisFile > len(validFiles)):
                print "Please enter an integer between 0 and", len(validFiles)
                continue
            elif runThisFile == len(validFiles):
                print "Exiting Program..."
                exitFlag = 1
                break
            break

        if exitFlag != 1:
            makeLammps = True
            plotGraph = True


            print "python CGer.py -x "+str(mode)+" -ml "+str(direc)+" -md "+str(validFiles[runThisFile][0])+" -mt "+str(validFiles[runThisFile][1])+' -d '+str(int(makeLammps))+' -g '+str(int(plotGraph))
            subprocess.call("python CGer.py -x "+str(mode)+" -ml "+str(direc)+" -md "+str(validFiles[runThisFile][0])+" -mt "+str(validFiles[runThisFile][1])+' -d '+str(int(makeLammps))+' -g '+str(int(plotGraph)), shell=True)
            continue
        else:
            break

    print "Battle control terminated."
            









    ### NOTE THAT THESE VALUES ARE CALIBRATED FOR 20 TOTAL CHAINS!!!
    #########################################################################################################
    ##### ---------==========  FIGURES FOR P3HT SIMILAR TO PLACES YOU CAN BUY ONLINE =========--------- #####
    #########################################################################################################
    ## Going for Mw = 54kDa, PDI ~ 1.7ish
    ## For 100 chains, totalMonomersPerChain = 280, testPolyDispIndex = 1.9
    ## For 40 chains: totalMonomersPerChain = 280, testPolyDispIndex = 1.9
    ## For 20 chains, the weight ends up being lower so compensate: totalMonomersPerChain = 300, testPolyDispIndex = 1.9. Note that 20 isn't enough for a proper statistical sample, so recommend to set checkBeforeMakingLammps = True and make sure the values of Mw and PDI are acceptable.
    ##
    ##
    ##
    ##
    ## polymerSource.com (3,000 < Mn < 31,200; 1.18 < PDisp < 1.45; https://secure.netsolhost.com/polymersource.com/productSearch.php?ID=P13182-3HT)
    ## This is also the same as the low-weight SigmaAldrich P3HT (http://www.sigmaaldrich.com/catalog/product/aldrich/698989?lang=en&region=GB)
    ## testPolyDispIndex = 1
    ## totalMonomersPerChain = 128
    #########################################################################################################
    ## solarisChem.com (50,000 < Mw < 70,000; 1.4 < PDisp < 1.6; http://www.solarischem.com/P3HT.html)
    ## testPolyDispIndex = 1.5
    ## totalMonomersPerChain = 384
    #########################################################################################################
    ## ossila.com (31,300 < Mw < 94,100; 1.7 < PDisp < 2.3; http://www.ossila.com/oled_opv_ofet_catalogue3/PCDTBT_P3HT_PCBM_PEDOT_PSS_for_organic_photovoltaics/M101-P3HT.php)
    ## This is also the same as the high-weight SigmaAldrich P3HT (http://www.sigmaaldrich.com/catalog/product/aldrich/698997?lang=en&region=GB)
    ## This is also the same as RiekeMetals Huang used in the Schwarz paper (http://www.riekemetals.com/component/content/article/81)
    ## testPolyDispIndex = 2.5
    ## totalMonomersPerChain = 320
    #########################################################################################################























    #         else: #submit == 'n'
    #             coresToUse = 0
    #             queue = None

    #         print "\n"
    #         t0 = T.time()
    #         subprocess.call('python analyseOrca.py -m '+str(morphologies[runThisFile])+' -c '+str(coresToUse)+' -q '+str(queue), shell=True)
    #         # subprocess.call('python analyseOrca.py -m '+str(morphologies[runThisFile]), shell=True)
    #         t1 = T.time()
    #         elapsedTime = float(t1) - float(t0)
    #         if elapsedTime < 60:
    #             timeunits = 'seconds.'
    #         elif elapsedTime < 3600:
    #             elapsedTime /= 60.0
    #             timeunits = 'minutes.'
    #         elif elapsedTime < 86400:
    #             elapsedTime /= 3600.0
    #             timeunits = 'hours.'
    #         else:
    #             elapsedTime /= 86400.0
    #             timeunits = 'days.'
    #         print "----------====================----------"
    #         print "analyseOrca.py calculations complete in %.1f %s." % (float(elapsedTime), str(timeunits))
    #         print "----------====================----------"

    # print "Battle control terminated."
