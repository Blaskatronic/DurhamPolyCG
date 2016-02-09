import os
import numpy as np
import pylab as P
import random as R
import csv
import subprocess
import argparse

def getDatList(direc):
    fileList = os.listdir(direc)
#    files = []
#    for element in fileList:
#        files.append(element[7:])
    datFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.dat'):
            datFiles.append(str(fileName))
    
    # for folder in dataFolders:
    #     path = str(direc)+str(folder)+'/conv/'
    #     convList = os.listdir(path)
    #     for iterationFile in convList:
    #         dataFiles.append(str(path)+str(iterationFile))
    # dataFiles = sorted(dataFiles)
    return datFiles


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


def generateMixture(chainLengthsDict, mass1, prop1, mass2, prop2):
    newMixture = []
    sample1 = chainLengthsDict[mass1]
    sample2 = chainLengthsDict[mass2]
    numberOfM1Chains = np.round(prop1/len(sample1))
    numberOfM2Chains = np.round(prop2/len(sample2))

    for i in range(numberOfM1Chains):
        chainLength = R.choice(sample1)
        sample1.remove(chainLength)
        newMixture.append(chainLength)

    for i in range(numberOfM2Chain):
        chainLength = R.choice(sample2)
        sample2.remove(chainLength)
        newMixture.append(chainLength)

    return newMixture   




def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    return locations
    



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m1", "--mass1", default=None, help="The weight averaged molecular weight (Mw) of the chains in the first sample to be mixed")
    parser.add_argument("-p1", "--prop1", default=None, help="The proportion of the final chain distribution that contains chains from the first sample")
    parser.add_argument("-m2", "--mass2", default=None, help="The weight averaged molecular weight (Mw) of the chains in the second sample to be mixed")
    parser.add_argument("-p2", "--prop2", default=None, help="The proportion of the final chain distribution that contains chains from the second sample")
    args = parser.parse_args()


    # TO DO: Add some foolproofing
    # If only one sample identified, then result is 100% mix of that. If -p1 given, then warning displayed.
    # If two samples are identified but no proportions, result is 50:50 mix. Warning displayed.
    # If either proportion is over 100%, result is a 50:50 mix (if 2 samples, 100% mix else) and warning is displayed.
    # If sum of proportions is over 100%, correct amount is added to/deducted from the second proportion. Warning displayed
    # If only one proportion is specified, then the rest is made up by the second sample. Notification displayed

    mass1 = int(args.mass1)
    mass2 = int(args.mass2)
    prop1 = float(args.prop1)
    prop2 = float(args.prop2)


    # if (mass1 != None) and (mass2 == None):
    #     print "Resultant sample will be 100% mix of sample with Mw =", args.mass1
    # elif (mass2 != None) and (mass1 == None):
    #     print "Resultant sample will be 100% mix of sample with Mw =", args.mass2
    # elif (mass1 == None) and (mass2 == None):
    #     raise SystemError('No arguements specified to make mixture')
    # else:
    #     if (prop1 > 1):



    direc = './inputDats/'
    validDats = getDatList(direc)
    validDats = sorted(validDats)
    print "Valid .dat files:", validDats


    chainLengthsDict = {}

    for datName in validDats:
        Mw = treatDatName(datName)
        masses, atoms, bonds, angles, dihedrals, impropers = loadDat(direc+datName)
        chainLengths = getChainLengths(atoms)
        print "Chain Lengths =", chainLengths
        chainLengthsDict[int(Mw)] = chainLengths


    customMixtureOfLengths = generateMixture(chainLengthsDict, mass1, prop1, mass2, prop2)

    print "Custom Mixture of", prop1, "Mw_"+str(mass1), "and", prop2, "Mw_"+str(mass2), ":"
    print customMixtureOfLengths
