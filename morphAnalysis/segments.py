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
 r  for datNumber in range(len(datFiles)):
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

def loadTrajectory(filename):
    timestepData = []
    timestepNos = []
    simVolData = []
    print "Reading in trajectory file:", filename
    originalFile = open(str(filename), 'r')
    trajData = originalFile.readlines()
    originalFile.close()

    nextLineTotalAtoms = 0
    nextLineAtomData = 0
    nextLineTimestep = 0
    nextLineSimVol = 0
    simVolLines = 0
    atomsRecorded = 0
    totalAtoms = 0
    lineNumber = 0

    for line in trajData:
        lineNumber += 1
        if nextLineSimVol == 1:
            if simVolLines == 3:
                simVolLines = 0
                nextLineSimVol = 0
            else:
                simVolData[-1].append([])
                simVolLines += 1
                dimensions = line.split(' ')
                for dim in dimensions:
                    simVolData[-1][-1].append(float(dim))
        if nextLineTimestep == 1:
            timestepNos.append(int(line[:-1]))
            nextLineTimestep = 0
        if nextLineTotalAtoms == 1:
            totalAtoms = int(line[:-1])
            nextLineTotalAtoms = 0
        if nextLineAtomData == 1:
            tempData = line[:-2].split(' ')
            for i in range(len(tempData)):
                # Unwrap periodicity and set to float/int as required
                # if (int(tempData[6]) != 0) or (int(tempData[7]) != 0) or (int(tempData[8]) != 0):
                #     print "--== BEFORE ==--"
                #     print "tempData =", tempData
                #     print "dim =", dim
                if (i == 3):
                    tempData[i] = float(tempData[i])+(int(tempData[i+3])*(2*float(dim)))
                elif (i == 4):
                    tempData[i] = float(tempData[i])+(int(tempData[i+3])*(2*float(dim)))
                elif (i == 5):
                    tempData[i] = float(tempData[i])+(int(tempData[i+3])*(2*float(dim)))
                else:
                    try:
                        tempData[i] = int(tempData[i])
                    except:
                        print "Line Number", lineNumber
                        print tempData
                        raise SystemError('FAIL')
                # if (int(tempData[6]) != 0) or (int(tempData[7]) != 0) or (int(tempData[8]) != 0):
                #     print "--== AFTER ==--"
                #     print "tempData =", tempData
                #     print "dim =", dim
            timestepData[-1].append(tempData)




            atomsRecorded += 1
        if 'ITEM: TIMESTEP' in line:
            nextLineTimestep = 1
        if 'ITEM: NUMBER OF ATOMS' in line:
            nextLineTotalAtoms = 1
        if 'ITEM: ATOMS' in line:
            timestepData.append([])
            nextLineAtomData = 1
        if 'ITEM: BOX BOUNDS' in line:
            simVolData.append([])
            nextLineSimVol = 1
        if atomsRecorded == totalAtoms:
            nextLineAtomData = 0
            atomsRecorded = 0

#    print timestepData[0]
    print "Timesteps = ", len(timestepData)
    print "Atoms =", len(timestepData[0])

    return timestepData, timestepNos, simVolData




def trimTrajectory(filename, molNumber):
    timestepData = []
    timestepNos = []
    simVolData = []
    print "Reading in trajectory file:", filename
    originalFile = open(str(filename), 'r')
    trajData = originalFile.readlines()
    originalFile.close()

    nextLineTotalAtoms = 0
    nextLineAtomData = 0
    nextLineTimestep = 0
    nextLineSimVol = 0
    simVolLines = 0
    atomsRecorded = 0
    totalAtoms = 0
    lineNumber = 0

    popList = []

    for line in trajData:
        lineNumber += 1
        if nextLineSimVol == 1:
            if simVolLines == 3:
                simVolLines = 0
                nextLineSimVol = 0
            else:
                simVolData[-1].append([])
                simVolLines += 1
                dimensions = line.split(' ')
                for dim in dimensions:
                    simVolData[-1][-1].append(float(dim))
        if nextLineTimestep == 1:
            timestepNos.append(int(line[:-1]))
            nextLineTimestep = 0
        if nextLineTotalAtoms == 1:
            totalAtoms = int(line[:-1])
            nextLineTotalAtoms = 0
        if nextLineAtomData == 1:
            tempData = line[:-2].split(' ')
            if tempData[2] != str(molNumber):
                popList.append(lineNumber-1)
            atomsRecorded += 1
        if 'ITEM: TIMESTEP' in line:
            nextLineTimestep = 1
        if 'ITEM: NUMBER OF ATOMS' in line:
            nextLineTotalAtoms = 1
        if 'ITEM: ATOMS' in line:
            timestepData.append([])
            nextLineAtomData = 1
        if 'ITEM: BOX BOUNDS' in line:
            simVolData.append([])
            nextLineSimVol = 1
        if atomsRecorded == totalAtoms:
            nextLineAtomData = 0
            atomsRecorded = 0


    popList.sort(reverse=True)

    for popValue in popList:
        try:
            trajData.pop(popValue)
        except:
            print popValue
            print len(trajData)
            raise SystemError('Cannot pop')


    singleMolTraj = open('./singleMolTraj.lammpstrj', 'w+')
    singleMolTraj.writelines(trajData)
    singleMolTraj.close()

    return

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

def sortByMolecule(atomsList):
    atomsByMolecule = []
    sortedAtomsList = sorted(atomsList, key=lambda x: x[2])
    totalMolecules = sortedAtomsList[-1][2]
    for i in range(totalMolecules):
        atomsByMolecule.append([])
    for atom in atomsList:
        atomsByMolecule[atom[2]-1].append(atom)
    return atomsByMolecule



def getAtomTypes(trajectoryData, mode):
    if mode == 'thio':
        atomType = 1
    elif mode == 'alk1':
        atomType = 2
    elif mode == 'alk2':
        atomType = 3
    atomsList = []
    for atom in trajectoryData:
        if atom[1] == atomType:
            atomsList.append(atom)
    sortedAtomsList = sortByMolecule(atomsList)
    return sortedAtomsList



def calculateSeparation(atom1, atom2):
    xdif = atom1[0] - atom2[0]
    ydif = atom1[1] - atom2[1]
    zdif = atom1[2] - atom2[2]
    return np.sqrt(xdif**2 + ydif**2 + zdif**2)

def normaliseVec(vector):
    return vector/float(np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2))

def piStackingDirection(angle):
    # Angle comes in as a 4-vector of:
    # [AtomNo, Type, Molecule, x, y, z, ix, iy, iz]
    # Normal direction to a plane defined by points P1, P2 and P3 =
    # (P3 - P1) vectorProduct (P2 - P1)
    # print "Atom =", atom
    vec1 = np.array([(angle[2][3] - angle[0][3]), (angle[2][4] - angle[0][4]), (angle[2][5] - angle[0][5])])
    vec2 = np.array([(angle[1][3] - angle[0][3]), (angle[1][4] - angle[0][4]), (angle[1][5] - angle[0][5])])
    normalVec = np.cross(vec1, vec2)
    normalVec = normaliseVec(normalVec)
    return normalVec


def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    return locations
    
def writeCSV(dataX, dataY, name):
    '''Appends a CSV file with X and Y Data'''
    filename = './'+name+'.csv'
    document = csv.writer(open(filename, 'a+'), delimiter = ',')
    document.writerow([dataX, dataY])


def plotMolecule(trajectoryData, segments, fileName):
    x = []
    y = []
    z = []

    colors = ['r', 'g', 'b', 'k']
    for segment in segments:
        x.append([])
        y.append([])
        z.append([])
        for atom in segment.iteritems():
            x[-1].append(atom[1][0])
            y[-1].append(atom[1][1])
            z[-1].append(atom[1][2])

    fig = P.figure()
    ax = p3.Axes3D(fig)
    for i in range(len(x)):
        colour = colors[i%4]
        for j in range(len(x[i])):
            ax.scatter(x[i][j], y[i][j], z[i][j], s = 20, c = colour)


    # colors = ['r', 'g', 'b', 'k']
    # for atom in trajectoryData:
    #     x.append(atom[3])#+(atom[6]*(abs(simVolData[0][0]) + abs(simVolData[0][1]))))
    #     y.append(atom[4])#+(atom[7]*(abs(simVolData[1][0]) + abs(simVolData[1][1]))))
    #     z.append(atom[5])#+(atom[8]*(abs(simVolData[2][0]) + abs(simVolData[2][1]))))

    # fig = P.figure()
    # ax = p3.Axes3D(fig)
    # for i in range(len(x)):
    #     ax.scatter(x[i], y[i], z[i], s = 20, c = colors[i%4])

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    P.savefig('./'+str(fileName)+'.png')
    # print "File saved as ./test"+filenameAddon+".png"




class molecule:
    def __init__(self, thiosInThisMolecule, alk1sInThisMolecule, alk2sInThisMolecule, bonds, toleranceAngle):
        self.allThios = self.makeDict(thiosInThisMolecule)
#        self.thiosInMorphology = self.makeDict(thiosInMorphology)
        self.thios = copy.deepcopy(self.allThios)
        self.alk1s = self.makeDict(alk1sInThisMolecule)
        self.alk2s = self.makeDict(alk2sInThisMolecule)
        self.thioThioBonds, self.thioAlk1Bonds, self.alk1Alk2Bonds = self.findRelevantBonds(bonds)
#        self.alkThioThioAngles = self.findRelevantAngles(angles)
        self.tolerance = toleranceAngle

        # Find Segments:
        self.segments = self.findSegments()


    def makeDict(self, atomsList):
        atomDict = {}
        for atom in atomsList:
            atomDict[atom[0]] = [atom[3], atom[4], atom[5]]
        return atomDict

    def findRelevantBonds(self, bonds):
        # Return bonds that are useful for this molecule. thioThioBonds are the backbone bonds, thioAlkBonds are the bonds between the thiophene and the first half of the alkyl chain (no longer used)
        thioThioBonds = []
        thioAlkBonds = []
        bondedAlk1s = []
        alk1Alk2Bonds = []
        for bond in bonds:
            if (bond[1] == 1): # Thio-thio bonds only
                if (bond[2] in self.thios) or (bond[3] in self.thios):
                    thioThioBonds.append([bond[2], bond[3]])
            elif (bond[1] == 2): # Thio-Alk1 bonds only
                if (bond[2] in self.thios):
                    thioAlkBonds.append([bond[2], bond[3]])
                    bondedAlk1s.append(bond[3])
                elif (bond[3] in self.thios):
                    thioAlkBonds.append([bond[3], bond[2]])
                    bondedAlk1s.append(bond[2])
        for bond in bonds:
            if (bond[1] == 3): # Alk1-Alk2 bonds only
                if (bond[2] in bondedAlk1s):
                    alk1Alk2Bonds.append([bond[2], bond[3]])
                elif (bond[3] in bondedAlk1s):
                    alk1Alk2Bonds.append([bond[3], bond[2]])
                
        return thioThioBonds, thioAlkBonds, alk1Alk2Bonds

    def findRelevantAngles(self, angles):
        alkThioThioAngles = []
        for angle in angles:
            if (angle[1] == 4): # Alk-Thio-Thio angles only
                alkThioThioAngles.append([angle[2], angle[3], angle[4]])
        return alkThioThioAngles
        

    def findBondedNeighbours(self, atom, currentSegment, previousSegment):
        bondedNeighbours = []
        for bond in self.thioThioBonds:
            if (bond[0] == atom) and (bond[1] not in currentSegment) and (bond[1] not in previousSegment):
                bondedNeighbours.append(bond[1])
            elif (bond[1] == atom) and (bond[0] not in currentSegment) and (bond[0] not in previousSegment):
                bondedNeighbours.append(bond[0])
        return bondedNeighbours

    def findStartingAtom(self):
        startingAtoms = {}
        for bond in self.thioThioBonds:
            if bond[0] not in startingAtoms:
                startingAtoms[bond[0]] = 1
            else:
                del startingAtoms[bond[0]]
            if bond[1] not in startingAtoms:
                startingAtoms[bond[1]] = 1
            else:
                del startingAtoms[bond[1]]
        return startingAtoms.keys()

    def findSegments(self):
        # Start with an atom that bonds in only one direction (end of chain)
        segmentMaster = [[]]
        moleculeEnds = self.findStartingAtom()
        atomUnderConsideration = moleculeEnds[0]
        atomUnderConsiderationCoords = self.thios[atomUnderConsideration]
        segmentMaster[-1].append(atomUnderConsideration)
        del self.thios[atomUnderConsideration]
        # The first run is a bit odd because we only have one bond and 2 atoms, so we have to recalculate this both before and then in the loop
        neighbouringThio = self.findBondedNeighbours(atomUnderConsideration, segmentMaster[-1], [0]) # Should only ever return one neighbouring thio as the "next one" along the chain

        firstAtomInSegment = True

        while len(self.thios) > 0:

            # print "---=== THIOS ===---"

            # for thio in self.thios.iteritems():
            #     print thio

            # print "---=============---\n"

            if len(segmentMaster) == 1:
                previousSegment = [0]
            else:
                previousSegment = segmentMaster[-2]

            neighbouringThio = self.findBondedNeighbours(atomUnderConsideration, segmentMaster[-1], previousSegment) # Should only ever return one neighbouring thio as the "next one" along the chain
            # print neighbouringThio
            if (firstAtomInSegment == True):
                axisVector = findAxis(atomUnderConsiderationCoords, self.thios[neighbouringThio[0]]) # Initial backbone axis
            separationVector = findAxis(atomUnderConsiderationCoords, self.thios[neighbouringThio[0]]) # Separation vector to next thiophene
            # print "Atom1 =", atomUnderConsiderationCoords, "...Atom2 =", self.thios[neighbouringThio[0]]
            # print "axisVector (should be the separation vector from previous iteration) =", axisVector
            # print "separationVector =", separationVector
            # crossProductMagnitude = np.linalg.norm(np.cross(axisVector, separationVector))
            # print "crossProductMagnitude =", crossProductMagnitude
            # separationAngle = np.arcsin(crossProductMagnitude)
            dotProduct = np.dot(axisVector, separationVector)
            # print "dotProduct =", dotProduct
            if abs(dotProduct-1.0) <= 1E-8:
                separationAngle = 0
            else:
                separationAngle = np.arccos(dotProduct)
            # if (firstAtomInSegment == False):
            #     csvFile = open('./separationAngles.csv', 'a+')
            #     csvWriter = csv.writer(csvFile, delimiter = ',')
            #     csvWriter.writerow([separationAngle])
            #     csvFile.close()
            # print "Separation Angle =", separationAngle
            if firstAtomInSegment == True:
                firstAtomInSegment = False
            else:
                axisVector = findAxis(atomUnderConsiderationCoords, self.thios[neighbouringThio[0]]) # Previous backbone vector, to be used in next loop iteration

            if abs(separationAngle) <= self.tolerance: # Tolerance angle
                pass
#                print "Adding", neighbouringThio[0], "to this segment. Current segment =", segmentMaster[-1]
            else:
                # print "Segment complete (angle =", str(separationAngle), "> 0.087), beginning new segment..."
                segmentMaster.append([])
                firstAtomInSegment = True
            segmentMaster[-1].append(neighbouringThio[0])

            atomUnderConsideration = neighbouringThio[0]
            atomUnderConsiderationCoords = self.thios[neighbouringThio[0]]
            del self.thios[neighbouringThio[0]]
            # print "\n"

        # print segmentMaster
        print "All thiophenes treated."
        # print "This molecule has", len(segmentMaster), "segments and looks like this:"
        # print segmentMaster

        return segmentMaster


    def findOrientations(self):
        # This is going to give me the normal vectors to the plane describing the alk1-atom-thio
        thioPlaneOrientationMaster = {} # The normal vector to the plane of the thiophene ring "AtomNo"
        alkylSideChainOrientationMaster = {} # The unit vector from the thiophene "AtomNo" to it's bonded Alk1
        neighbouringThiosCoords = {} # The position(s) of the neighbouring thiophene(s) bonded to the thiophene "AtomNo"
        bondedAlk1sCoords = {} # The position of the Alk1 bonded to the thiophene "AtomNo"
        bondedAlk2sCoords = {} # The position of the Alk2 bonded to the thiophene "AtomNo"
        for segment in self.segments:
            previousBond = 0
            for atomNo in segment:
                thioCoords = self.allThios[atomNo]
                # First, find the separation vector between the two neighbouring thios. This will match the plane of the thiophene ring in question
                # If only one neighbour, match the plane to the separation vector between the thio in question and it's neighbour
                bondedThiophenes = []
                for bond in self.thioThioBonds:
                    if atomNo == bond[0]:
                        bondedThiophenes.append(bond[1])
                    elif atomNo == bond[1]:
                        bondedThiophenes.append(bond[0])

                if len(bondedThiophenes) == 2:
                    firstThioCoords = self.allThios[bondedThiophenes[0]]
                    secondThioCoords = self.allThios[bondedThiophenes[1]]
                    # NB THESE HAVE TO BE DEEPCOPIED TO PREVENT A RECURSION ERROR
                    neighbouringThiosCoords[atomNo] = [copy.deepcopy(firstThioCoords), copy.deepcopy(secondThioCoords)]
                elif len(bondedThiophenes) == 1:
                    firstThioCoords = self.allThios[atomNo]
                    secondThioCoords = self.allThios[bondedThiophenes[0]]
                    # NB THESE HAVE TO BE DEEPCOPIED TO PREVENT A RECURSION ERROR
                    neighbouringThiosCoords[atomNo] = [copy.deepcopy(secondThioCoords)]
                else:
                    print "Current Atom Number =", atomNo
                    print "Bonded Thiophenes =", bondedThiophenes
                    raise SystemError('Bond Calculation Error - more than 2 or zero bonded thiophenes')

                adjThiosSeparationVector = [secondThioCoords[0] - firstThioCoords[0], secondThioCoords[1] - firstThioCoords[1], secondThioCoords[2] - firstThioCoords[2]]

                adjThiosSeparationVector = normaliseVec(np.array(adjThiosSeparationVector))


                # foundThioThioBond = False
                # for bond in self.thioThioBonds:
                #     if (bond != previousBond) and (atomNo == bond[0]):
                #         previousBond = bond
                #         foundThioThioBond = True
                #         adjThioCoords = self.allThios[bond[1]]
                #         break
                #     elif (bond != previousBond) and (atomNo == bond[1]):
                #         previousBond = bond
                #         adjThioCoords = self.allThios[bond[0]]
                #         foundThioThioBond = True
                #         break
                # if foundThioThioBond == False:
                #     bond = previousBond
                #     if (atomNo == bond[0]):
                #         foundThioThioBond = True
                #         adjThioCoords = self.allThios[bond[1]]
                #     elif (atomNo == bond[1]):
                #         foundThioThioBond = True
                #         adjThioCoords = self.allThios[bond[0]]


                # Now get the Thio-Alk1 bond

                for bond in self.thioAlk1Bonds:
                    if atomNo == bond[0]:
                        alk1Coords = self.alk1s[bond[1]]
                        alk1AtomNo = bond[1]
                        break
                    elif atomNo == bond[1]:
                        alk1Coords = self.alk1s[bond[0]]
                        alk1AtomNo = bond[0]
                        break

                for bond in self.alk1Alk2Bonds:
                    if alk1AtomNo == bond[0]:
                        alk2Coords = self.alk2s[bond[1]]
                        break
                    elif alk1AtomNo == bond[1]:
                        alk2Coords = self.alk2s[bond[0]]
                        break

                bondedAlk1sCoords[atomNo] = alk1Coords
                bondedAlk2sCoords[atomNo] = alk2Coords

                alk1SeparationAxis = findAxis(thioCoords, alk1Coords)

                # The normal to the thio ring is the cross product of the alk1SeparationAxis and the adjThiosSeparationVector
                normalToThioRing = np.cross(alk1SeparationAxis, adjThiosSeparationVector)
                normalToThioRing = normaliseVec(normalToThioRing)

                # thioSeparationAxis = findAxis(thioCoords, adjThioCoords)
                # normalToPlane = np.cross(alk1SeparationAxis, thioSeparationAxis)
#                normalToPlaneMagnitude = np.sqrt(((normalToPlane[0])**2) + ((normalToPlane[1])**2) + ((normalToPlane[2])**2))
 #               normalToPlane /= normalToPlaneMagnitude
                # normalToPlane = normaliseVec(normalToPlane)
                # thioNormalOrientationMaster[atomNo] = [normalToPlane, adjThioCoords]
                thioPlaneOrientationMaster[atomNo] = normalToThioRing
                alkylSideChainOrientationMaster[atomNo] = alk1SeparationAxis
        return thioPlaneOrientationMaster, alkylSideChainOrientationMaster, neighbouringThiosCoords, bondedAlk1sCoords, bondedAlk2sCoords



    def returnSegments(self):
        segmentReturn = []
        for segment in self.segments:
            segmentReturn.append({})
            for atom in segment:
                segmentReturn[-1][atom] = self.allThios[atom]
        return segmentReturn

    def returnOrientations(self):
        thioOrientations, alkylSideChainOrientations, neighbouringThiosCoords, bondedAlk1sCoords, bondedAlk2sCoords = self.findOrientations()
        return thioOrientations, alkylSideChainOrientations, neighbouringThiosCoords, bondedAlk1sCoords, bondedAlk2sCoords


        
def getFilename(dataFile, moleculeNo):

    underscores = findIndex(dataFile, '_')

    moleculeNumber = str(moleculeNo+1)
    while len(moleculeNumber) < 3:
        moleculeNumber = '0'+moleculeNumber

    filename = dataFile[:underscores[1]]+dataFile[underscores[-2]:-4]+'_mol'+moleculeNumber

    return filename





# -------------==================== Getting Pairs of Segments  ======================-----------------

def calcSeparation(coords1, coords2):
    # Calculates the physical separation between two positions of the form
    # [x, y, z]
    separation = np.sqrt((coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 + (coords1[2]-coords2[2])**2)
    return separation


def findAxis(atom1, atom2):
    xSep = atom2[0] - atom1[0]
    ySep = atom2[1] - atom1[1]
    zSep = atom2[2] - atom1[2]
    axisVector = normaliseVec(np.array([xSep, ySep, zSep]))
    return axisVector


def findPeriodicSegmentsList(segmentsMaster, maximumHoppingDistance, simVolData):
    periodicMaster = []
    rollingSegmentNumber = 0
    atomsInImages = 0
    for segment in segmentsMaster:
        rollingSegmentNumber += 1
        requiresImageSegment = False
        imageCoords = [0, 0, 0]
        for atom in segment.iteritems():
            # Like in LAMMPS, generate an "image position" list with [ix, iy, iz] being -1 or 1 depending on whether it is within maximumHoppingDistance of a simulation boundary
            # xCoords, yCoords and zCoords should technically check Alk1 and Alk2 too!
            for xCoords in [atom[1][0], atom[1][8][0], atom[1][9][0]]: # xCoords of thiophene, alk1 and alk2
                if (xCoords <= (simVolData[0][0]+maximumHoppingDistance)):
                    requiresImageSegment = True
                    if imageCoords[0] != 1:
                        imageCoords[0] = 1
                elif (xCoords >= (simVolData[0][1]-maximumHoppingDistance)):
                    requiresImageSegment = True
                    if imageCoords[0] != -1:
                        imageCoords[0] = -1

            for yCoords in [atom[1][1], atom[1][8][1], atom[1][9][1]]: # yCoords of thiophene, alk1 and alk2
                if (yCoords <= (simVolData[1][0]+maximumHoppingDistance)):
                    requiresImageSegment = True
                    if imageCoords[1] != 1:
                        imageCoords[1] = 1
                elif (yCoords >= (simVolData[1][1]-maximumHoppingDistance)):
                    requiresImageSegment = True
                    if imageCoords[1] != -1:
                        imageCoords[1] = -1

            for zCoords in [atom[1][2], atom[1][8][2], atom[1][9][2]]: # zCoords of thiophene, alk1 and alk2
                if (zCoords <= (simVolData[2][0]+maximumHoppingDistance)):
                    requiresImageSegment = True
                    if imageCoords[2] != 1:
                        imageCoords[2] = 1
                elif (zCoords >= (simVolData[2][1]-maximumHoppingDistance)):
                    requiresImageSegment = True
                    if imageCoords[2] != -1:
                        imageCoords[2] = -1

        if requiresImageSegment == True:
            periodicSegment = {}
            imageSegment = copy.deepcopy(segment)
            xImageDistance = simVolData[0][1] - simVolData[0][0]
            yImageDistance = simVolData[1][1] - simVolData[1][0]
            zImageDistance = simVolData[2][1] - simVolData[2][0]
            for imageAtom in imageSegment.iteritems():
                # Treat ThioX
                imageAtom[1][0] += imageCoords[0]*xImageDistance
                # Treat ThioY
                imageAtom[1][1] += imageCoords[1]*yImageDistance
                # Treat ThioZ
                imageAtom[1][2] += imageCoords[2]*zImageDistance
                # Treat neighbouring thiophenes (there can be one or two of these)
                for neighbourThiophene in imageAtom[1][7]:
                    neighbourThiophene[0] += imageCoords[0]*xImageDistance
                    neighbourThiophene[1] += imageCoords[1]*yImageDistance
                    neighbourThiophene[2] += imageCoords[2]*zImageDistance
                # Treat Bonded Alk1
                imageAtom[1][8][0] += imageCoords[0]*xImageDistance
                imageAtom[1][8][1] += imageCoords[1]*yImageDistance
                imageAtom[1][8][2] += imageCoords[2]*zImageDistance
                # Treat Bonded Alk2
                imageAtom[1][9][0] += imageCoords[0]*xImageDistance
                imageAtom[1][9][1] += imageCoords[1]*yImageDistance
                imageAtom[1][9][2] += imageCoords[2]*zImageDistance
                periodicSegment[imageAtom[0]] = imageAtom[1]
            periodicMaster.append(periodicSegment)
    return periodicMaster


            




def findSegmentPairs(segmentsMaster, atomMaster, maximumHoppingDistance, simVolData):
    # len(segmentsMaster) = number of segments
    print simVolData
    periodicMaster = findPeriodicSegmentsList(segmentsMaster, maximumHoppingDistance, simVolData) # The same structure as segmentsMaster, but only including segments that are not in the initial bounded volume (i.e. are periodically linked to the main system)
    # These periodic segments have the same segment number as in the main, non-image volume, so as long as we use that name for the orca input files, but use the coordinates given in periodicMaster, then we should be ok to calculate transfer integrals.
    # Because of this, treat `periodic pairs' separately to the main segment pairs to ensure we don't get mixed up.
    periodicMasterDict = {} # Useful to sort the list by the segment number for the periodic segments because it no longer lines up with the list indices
    periodicAtomMaster = {}
    for periodicSegment in periodicMaster:
        for atom in periodicSegment.iteritems():
            periodicAtomMaster[atom[0]] = atom[1]

    segmentPairs = []
    periodicSegmentPairs = []
    # Just look at first one first
    for i in range(len(segmentsMaster)):
        segment = segmentsMaster[i]
        nearbySegments = []
        nearbyPeriodicSegments = []
        for currentAtom in segment.iteritems():
            currentAtomCoords = [currentAtom[1][0], currentAtom[1][1], currentAtom[1][2]]
            currentSegment = currentAtom[1][4]
            for compareAtom in atomMaster.iteritems():
                if (compareAtom[1][4] in nearbySegments) or (compareAtom[1][4] == currentSegment) or (currentAtom[0] == compareAtom[0]):
                    continue
                compareAtomCoords = [compareAtom[1][0], compareAtom[1][1], compareAtom[1][2]]
                separation = calcSeparation(currentAtomCoords, compareAtomCoords)
                if separation >= maximumHoppingDistance: # Arbitrarily discard atoms more than maximumHoppingDistance angstroms away
                    continue
                # If it's not continued already:
                # Atom's home segment is not classed as nearby to currentSegment and is within 10 angstroms
                nearbySegments.append(compareAtom[1][4])
            for compareAtom in periodicAtomMaster.iteritems():
                if (compareAtom[1][4] in nearbySegments) or (compareAtom[1][4] == currentSegment) or (currentAtom[0] == compareAtom[0]):
                    continue
                compareAtomCoords = [compareAtom[1][0], compareAtom[1][1], compareAtom[1][2]]
                separation = calcSeparation(currentAtomCoords, compareAtomCoords)
                if separation >= maximumHoppingDistance:
                    continue
                if compareAtom[1][4] not in nearbyPeriodicSegments:
                    nearbyPeriodicSegments.append(compareAtom[1][4])

        if [currentSegment] not in segmentPairs:
            segmentPairs.append([currentSegment])

        # Remove the current segment from the nearby list
        for i in range(len(nearbySegments)):
            if nearbySegments[i] == currentSegment:
                continue
            if currentSegment < nearbySegments[i]:
                segmentPair = [currentSegment, nearbySegments[i]]
            else:
                segmentPair = [nearbySegments[i], currentSegment]
            if segmentPair not in segmentPairs:
                segmentPairs.append(segmentPair)

        for i in range(len(nearbyPeriodicSegments)):
            if nearbyPeriodicSegments[i] == currentSegment:
                continue
            periodicSegmentPair = [currentSegment, nearbyPeriodicSegments[i]] # Always real segment first, periodic second
            if periodicSegmentPair not in periodicSegmentPairs:
                periodicSegmentPairs.append(periodicSegmentPair)

        print "Current Segment =", currentSegment, "Nearby `real' Segments =", nearbySegments, "Nearby `periodic' Segments =", nearbyPeriodicSegments

    # print "SegmentPairs =", segmentPairs
    print "\nLen SegmentPairs =", len(segmentPairs)
    print "Len periodicSegmentPairs =", len(periodicSegmentPairs)
    print "Job done."

    for segment in periodicMaster:
        segmentNumber = segment.items()[0][1][4] # The number for this segment
        periodicMasterDict[segmentNumber] = segment

    return segmentPairs, periodicSegmentPairs, periodicMasterDict



# -------------==================== ======================-----------------




# -------------==================== Creating ORCA Inputs  ======================-----------------

class ORCAInput:
    def __init__(self, segmentsMaster, segmentPair, datName, rotateThio, rotateAlk1, rotateAlk2, justPlotThios, periodicSegmentPair = False, periodicMasterDict = None):
        self.rotateThio = rotateThio
        self.rotateAlk1 = rotateAlk1
        self.rotateAlk2 = rotateAlk2
        self.justPlotThios = justPlotThios
        self.segmentPair = segmentPair
        self.datName = datName[:-4]
        self.segmentAtoms = []
        if periodicSegmentPair == False:
            self.segmentAtoms.append(segmentsMaster[segmentPair[0]-1])
            COMSeg1 = self.determineCOMFromCG(segmentsMaster[segmentPair[0]-1])
            COMSeg2 = None
            if len(segmentPair) == 2:
                COMSeg2 = self.determineCOMFromCG(segmentsMaster[segmentPair[1]-1])
                self.segmentAtoms.append(segmentsMaster[segmentPair[1]-1])
            self.writeCOMCSV(COMSeg1, COMSeg2)
        else:
            self.segmentAtoms.append(segmentsMaster[segmentPair[0]-1])
            self.segmentAtoms.append(periodicMasterDict[segmentPair[1]])
            COMSeg1 = self.determineCOMFromCG(segmentsMaster[segmentPair[0]-1])
            COMSeg2 = self.determineCOMFromCG(periodicMasterDict[segmentPair[1]])
            self.writeCOMCSV(COMSeg1, COMSeg2)
        self.monomerAtoms = self.read3HTTemplate() # These are the coordinates of the atoms representing a CG site at [0, 0, 0]

    def writeCOMCSV(self, COMData1, COMData2):
        # Read the current data into a dictionary to make sure that we're not re-writing the same data
        currentSegData = {}
        if COMData2 == None:
            COMCSVName = './orcaOutputs/'+self.datName+'/SegmentSingleCOM.csv'
            csvReadOnly = open(COMCSVName, 'a+')
            csvReadOnly.seek(0)
            csvData = csv.reader(csvReadOnly, delimiter = ',')
            for row in csvData:
                currentSegData[row[0]] = 1
            csvReadOnly.close()
            currentSegment = COMData1[0]
            if (currentSegment not in currentSegData):
                csvFile = open(COMCSVName, 'a')
                csvWriter = csv.writer(csvFile, delimiter = ',')
                csvWriter.writerow([COMData1[0], COMData1[1][0], COMData1[1][1], COMData1[1][2]])
                csvFile.close()
        else:
            COMCSVName = './orcaOutputs/'+self.datName+'/SegmentPairCOM.csv'
            csvReadOnly = open(COMCSVName, 'a+')
            csvReadOnly.seek(0)
            csvData = csv.reader(csvReadOnly, delimiter = ',')
            for row in csvData:
                currentSegData[(row[0], row[1])] = 1
            csvReadOnly.close()
            currentSegmentPair = (COMData1[0], COMData2[0])
            if (currentSegmentPair not in currentSegData):
                separation = calcSeparation([COMData1[1][0], COMData1[1][1], COMData1[1][2]], [COMData2[1][0], COMData2[1][1], COMData2[1][2]])
                csvFile = open(COMCSVName, 'a')
                csvWriter = csv.writer(csvFile, delimiter = ',')
                # Write the row as [seg1ID, seg2ID, seg1X, seg1Y, seg1Z, seg2X, seg2Y, seg2Z, separation]
                csvWriter.writerow([COMData1[0], COMData2[0], COMData1[1][0], COMData1[1][1], COMData1[1][2], COMData2[1][0], COMData2[1][1], COMData2[1][2], separation])
                csvFile.close()

    def determineCOMFromCG(self, inputAtoms):
        # SegmentsMaster-like object comes in
        massOfThio = 81.11657
        massOfAlk1 = 42.08127
        massOfAlk2 = 43.08924
        massWeightedX = 0.
        massWeightedY = 0.
        massWeightedZ = 0.
        totalMass = 0.
        for atomData in inputAtoms.values():
            currentSegment = atomData[4]

            massWeightedX += atomData[0]*massOfThio
            massWeightedY += atomData[1]*massOfThio
            massWeightedZ += atomData[2]*massOfThio
            totalMass += massOfThio

            massWeightedX += atomData[8][0]*massOfAlk1
            massWeightedY += atomData[8][1]*massOfAlk1
            massWeightedZ += atomData[8][2]*massOfAlk1
            totalMass += massOfAlk1

            massWeightedX += atomData[9][0]*massOfAlk2
            massWeightedY += atomData[9][1]*massOfAlk2
            massWeightedZ += atomData[9][2]*massOfAlk2
            totalMass += massOfAlk2
        return [currentSegment, np.array([massWeightedX/float(totalMass), massWeightedY/float(totalMass), massWeightedZ/float(totalMass)])]


    def determineCOM(self, inputAtoms, thioOnly = False, alk1Only = False, alk2Only = False):
        massWeightedX = 0.
        massWeightedY = 0.
        massWeightedZ = 0.
        totalMass = 0.
        atoms = []
        if thioOnly == True: # These are the lines that contain thiophene atoms in the .xyz
            atoms.append(inputAtoms[0]) # S1
            atoms.append(inputAtoms[1]) # C2
            atoms.append(inputAtoms[3]) # C4
            atoms.append(inputAtoms[4]) # C5
            atoms.append(inputAtoms[5]) # C6
            atoms.append(inputAtoms[7]) # H8
        elif alk1Only == True: # These are the lines that contain alk1 atoms in the .xyz
            atoms.append(inputAtoms[8]) # C9
            atoms.append(inputAtoms[10]) # C11
            atoms.append(inputAtoms[13]) # C14
            atoms.append(inputAtoms[12]) # H13
            atoms.append(inputAtoms[2]) # H3
            atoms.append(inputAtoms[7]) # H6
            atoms.append(inputAtoms[15]) # H16
            atoms.append(inputAtoms[9]) # H10
            atoms.append(inputAtoms[17]) # H18
        elif alk2Only == True: # These are the lines that contain alk2 atoms in the .xyz
            atoms.append(inputAtoms[16])
            atoms.append(inputAtoms[11])
            atoms.append(inputAtoms[20])
            atoms.append(inputAtoms[19])
            atoms.append(inputAtoms[14])
            atoms.append(inputAtoms[22])
            atoms.append(inputAtoms[21])
            atoms.append(inputAtoms[18])
            atoms.append(inputAtoms[23])
            atoms.append(inputAtoms[24])
        else:
            atoms = inputAtoms
        for atom in atoms:
            if atom[0] == 0: # Determining chain COM, therefore mass factor is irrelevent
                mass = 1.0                
            elif atom[0] == 'S':
                mass = 32.065
            elif atom[0] == 'C':
                mass = 12.0107
            elif atom[0] == 'H':
                mass = 1.00794

            massWeightedX += atom[1][0]*mass
            massWeightedY += atom[1][1]*mass
            massWeightedZ += atom[1][2]*mass
            totalMass += mass
        return np.array([massWeightedX/float(totalMass), massWeightedY/float(totalMass), massWeightedZ/float(totalMass)])

    def read3HTTemplate(self):
        monomerFile = open('./templates/3HT_topdown.xyz', 'r')
        monomer = monomerFile.readlines()
        monomerFile.close()
        # Check centrepoint
        atoms = []
        for line in monomer[2:-1]:
            tempLine = []
            coordinates = []
            for character in line.split(' '):
                if len(character) != 0:
                    if len(character) == 1: # This is the element of the atom
                        tempLine.append(str(character))
                    else: # This is a coordinate
                        try:
                            coordinates.append(float(character))
                        except: # "\n" is breaking it
                            coordinates.append(float(character[:-2]))
            tempLine.append(coordinates)
            atoms.append(tempLine)

        COM = self.determineCOM(atoms, thioOnly = True)
        for atom in atoms:
            coordArray = np.array(atom[1])
            coordArray -= COM
            atom[1] = list(coordArray)

        return atoms


    def createName(self):
        # Name was too complicated, made it easier
        if (len(segmentPair) == 1):
            segNumber = str(segmentPair[0])
            while len(segNumber) < 4:
                segNumber = '0'+segNumber
            self.name = 'seg'+segNumber+'.inp'
        elif (len(segmentPair) == 2):
            segNumber1 = str(segmentPair[0])
            while len(segNumber1) < 4:
                segNumber1 = '0'+segNumber1
            segNumber2 = str(segmentPair[1])
            while len(segNumber2) < 4:
                segNumber2 = '0'+segNumber2
            self.name = 'seg'+segNumber1+'_seg'+segNumber2+'.inp'
        else:
            print "Segment Pair =", segmentPair
            raise SystemError("Segment Pair length is either 0 or > 2")


        # # Want a name kinda like C1_LengthOfChain1_C2_LengthOfChain2_X_XSeparationFromCentrePoint_Y_YSep_Z_ZSep.inp
        # C1Length = len(self.segmentAtoms[0])
        # C2Length = len(self.segmentAtoms[1])
        # # Determine COMs of each chain
        # chain1 = []
        # chain2 = []
        # for atom in self.segmentAtoms[0].iteritems():
        #     chain1.append([0, [atom[1][0], atom[1][1], atom[1][2]]])
        # for atom in self.segmentAtoms[1].iteritems():
        #     chain2.append([0, [atom[1][0], atom[1][1], atom[1][2]]])
        # COM1 = self.determineCOM(chain1)
        # COM2 = self.determineCOM(chain2)
        # separation = COM2-COM1
        # self.name = 'C1_'+str(len(chain1))+'_C2_'+str(len(chain2))+'_X_%.1f_Y_%.1f_Z_%.1f' % (float(separation[0]), float(separation[1]), float(separation[2]))


    def makeDirs(self):
        # Make the correct directory
        dirs = os.listdir('./orcaInputs')
        if self.datName not in dirs:
            os.makedirs('./orcaInputs/'+self.datName)
        filesList = os.listdir('./orcaInputs/'+self.datName)
        self.fullPath = './orcaInputs/'+self.datName+'/'+self.name
        if self.name in filesList:
            return True
        else:
            return False

    def determineThioToAlk1Vector(self, atoms):
        #thioPosn = [0, 0, 0] because we moved the monomer to the middle when we read it
        thioPosn = self.determineCOM(atoms, thioOnly = True)
        alk1Posn = self.determineCOM(atoms, alk1Only = True)
        # print "ThioPosn =", thioPosn
        # print "alk1Posn =", alk1Posn
        alk1Vector = normaliseVec(alk1Posn-thioPosn)
        return alk1Vector



    def determineNormalToThiopheneRing(self, monomerAtoms):
        # Plane is determined by two vectors, one from the Sulphur (S1 - atom 1 in the .xyz) to a neighbouring Carbon (C4 - atom 4 in the .xyz),
        # and the other between the sulphur (S1 - atom 1) to a carbon on the opposite side of the ring (C5 - atom 5 in the .xyz)
#        sulphurAtom = monomerAtoms[0]

        firstCarbon = monomerAtoms[1]
        neighbourCarbon = monomerAtoms[4]
        oppositeCarbon = monomerAtoms[3]
        
        vec1 = findAxis(firstCarbon[1], neighbourCarbon[1])
        vec2 = findAxis(firstCarbon[1], oppositeCarbon[1])


        normalToRing = np.cross(vec2, vec1)
        # normalToRingMagnitude = np.sqrt(((normalToRing[0])**2) + ((normalToRing[1])**2) + ((normalToRing[2])**2))
        # normalToRing /= normalToRingMagnitude
        normalToRing = normaliseVec(normalToRing)

        return normalToRing


    # def determineNormalToThiopheneRing(self, monomerAtoms, thioAlk1Vector):
    #     # Plane is determined by two vectors, one from the thiophene-alk1 vector (given)
    #     # And the other between the atoms C2 and C4 (which are the carbons that bond to the adjacent monomers)
    #     print monomerAtoms[1], monomerAtoms[3]
    #     leftCarbon = monomerAtoms[1][1] # C2
    #     rightCarbon = monomerAtoms[3][1] # C4

    #     crossRingVector = findAxis(leftCarbon, rightCarbon)

    #     print crossRingVector, "\n"

    #     normalToRing = np.cross(crossRingVector, thioAlk1Vector)
    #     normalToRing = normaliseVec(normalToRing)

    #     return normalToRing

    def rotationMatrixAroundAxis(self, theta, axis):
        # Rotation matrix obtained from http://goo.gl/RkW80
        R = np.matrix([[ np.cos(theta)+(axis[0]**2)*(1-np.cos(theta)), axis[0]*axis[1]*(1-np.cos(theta))-axis[2]*np.sin(theta), axis[0]*axis[2]*(1-np.cos(theta))+axis[1]*np.sin(theta) ],
                       [ axis[1]*axis[0]*(1-np.cos(theta))+axis[2]*np.sin(theta), np.cos(theta)+(axis[1]**2)*(1-np.cos(theta)), axis[1]*axis[2]*(1-np.cos(theta))-axis[0]*np.sin(theta) ],
                       [ axis[2]*axis[0]*(1-np.cos(theta))-axis[1]*np.sin(theta), axis[2]*axis[1]*(1-np.cos(theta))+axis[0]*np.sin(theta), np.cos(theta)+(axis[2]**2)*(1-np.cos(theta)) ]])
        return R


    def executeRotation(self, testRotationAngle, rotationAxis, leftCarbon, rightCarbon):
        leftCarbon = copy.deepcopy(leftCarbon)
        rightCarbon = copy.deepcopy(rightCarbon)
        rotationMatrix = self.rotationMatrixAroundAxis(testRotationAngle, rotationAxis)
        newLeftCarbonCoords = np.array(np.transpose(rotationMatrix*np.transpose(np.matrix(leftCarbon))))[0].tolist()
        newRightCarbonCoords = np.array(np.transpose(rotationMatrix*np.transpose(np.matrix(rightCarbon))))[0].tolist()
        return newLeftCarbonCoords, newRightCarbonCoords


    def calcSeparationAngle(self, vec1, vec2):
        vec1 = normaliseVec(vec1)
        vec2 = normaliseVec(vec2)
        dotProduct = np.dot(vec1, vec2)
        theta = np.arccos(dotProduct)
        theta = abs(theta)
        if theta > np.pi:
            theta = theta-np.pi
        return theta



    def thioInPlaneRotation(self, atoms, adjThioCoords, rotationAxis):
        # print "\n"
        totalAngleTurnedBy = 0
        testRotationAngle = np.pi/180.
        # Target vector is the vector between the two neighbour thiophenes, or just the vector to the adjacent thio if only one bonded neighbour
        if len(adjThioCoords) == 2:
            targetVector = findAxis(adjThioCoords[0], adjThioCoords[1])
        else:
            # Ends of chains don't really work because they only have one neighbour. As such, let's ignore them entirely and do not rotate
            return np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        # Current vector is the vector between the left and right carbons in the chain
        leftCarbon = copy.deepcopy(atoms[3][1])
        rightCarbon = copy.deepcopy(atoms[1][1])
        # print "Left Carbon", leftCarbon, "Right Carbon", rightCarbon
        currentVector = findAxis(leftCarbon, rightCarbon)
        currentTheta = self.calcSeparationAngle(currentVector, targetVector)
        # To minimise, perform rotation. If arccos of dot product > pi, subtract pi. (just means going opposite way so that's fine). Minimise theta.
        # First, perform a test rotation to see if we're rotating the right way
        # print "Target vector =", targetVector
        # print "Current vector =", currentVector
        # print "Beginning test rotation..."
        newLeftCarbonCoords, newRightCarbonCoords = self.executeRotation(testRotationAngle, rotationAxis, leftCarbon, rightCarbon)
        newVector = findAxis(newLeftCarbonCoords, newRightCarbonCoords)
        newTheta = self.calcSeparationAngle(newVector, targetVector)
        # print "Current Theta =", currentTheta, "new Theta =", newTheta
        if newTheta > currentTheta:
            # print "Rotating the wrong way! Reversing test rotation..."
            # Rotating the wrong way!
            testRotationAngle = -testRotationAngle
            newLeftCarbonCoords, newRightCarbonCoords = self.executeRotation(testRotationAngle, rotationAxis, leftCarbon, rightCarbon)
            newVector = findAxis(newLeftCarbonCoords, newRightCarbonCoords)
            newTheta = self.calcSeparationAngle(newVector, targetVector)
            if newTheta > currentTheta:
                print "Got an issue, rotated the opposite way but it didn't help."
                print "Likely means that we are already at the optimal position"
                return np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        # newTheta is now < currentTheta
        # print "---=== NEW VECTOR =", newVector, "===---"
        while newTheta <= currentTheta:
            # print "Beginning of loop: New Theta =", newTheta, "Current Theta =", currentTheta
            currentTheta = copy.deepcopy(newTheta)
            currentVector = copy.deepcopy(newVector)
            totalAngleTurnedBy += testRotationAngle
            # print "Executing another rotation (currentTheta =", str(currentTheta)+", totalAngleTurnedBy =", str(totalAngleTurnedBy)+")"
            newLeftCarbonCoords, newRightCarbonCoords = self.executeRotation(testRotationAngle, rotationAxis, leftCarbon, rightCarbon)
            # print "Left Carbon =", leftCarbon, "newLeftCarbon =", newLeftCarbonCoords
            # print "Right Carbon =", rightCarbon, "newRightCarbon =", newRightCarbonCoords
            newVector = findAxis(newLeftCarbonCoords, newRightCarbonCoords)
            newTheta = self.calcSeparationAngle(newVector, targetVector)
            # print "NEW VECTOR =", newVector
            # print "End of loop: Current theta is still", currentTheta, " but now newTheta =", newTheta, "(which has to be < currentTheta)"
            leftCarbon = newLeftCarbonCoords
            rightCarbon = newRightCarbonCoords

        rotMatrix = self.rotationMatrixAroundAxis(totalAngleTurnedBy, rotationAxis)


        # print "NewLeftCarbon", newLeftCarbonCoords, "NewRightCarbon", newRightCarbonCoords
        # print "NewVector =", newVector

        return rotMatrix
                


    def rotationMatrix(self, vector1, vector2):
        # A function to return the rotation matrix around the origin that maps vector1 to vector 2
        crossProduct = np.cross(vector1, vector2)
        sinAngle = np.sqrt(((crossProduct[0]**2) + ((crossProduct[1])**2) + ((crossProduct[2])**2)))
        cosAngle = np.dot(vector1, vector2)

        skewMatrix = np.matrix([[0, -crossProduct[2], crossProduct[1]], [crossProduct[2], 0, -crossProduct[0]], [-crossProduct[1], crossProduct[0], 0]])
        skewMatrixSquared = skewMatrix * skewMatrix
        
        rotMatrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix + skewMatrixSquared*((1 - cosAngle)/(sinAngle**2))

        return rotMatrix



    def manipulateMonomer(self, COMThio, CGThioPlaneNormal, CGAdjThioCoords, CGThioAlk1, COMAlk1, COMAlk2):
        alk1ManipulatedMonomer = []
        alk1AndThioManipulatedMonomer = []
        alk1AndThioZManipulatedMonomer = []
        thioAlk1CompletedMonomer = []
        thioAlk1Alk2CompletedMonomer = []
        monomerFinal = []

        relativeAlk1Posn = [(COMAlk1[0] - COMThio[0]), (COMAlk1[1] - COMThio[1]), (COMAlk1[2] - COMThio[2])] # Alk1 posn with thio as origin
        relativeAlk2Posn = [(COMAlk2[0] - COMThio[0]), (COMAlk2[1] - COMThio[1]), (COMAlk2[2] - COMThio[2])] # Alk1 posn with thio as origin

        templateAlk1Vector = self.determineThioToAlk1Vector(self.monomerAtoms)

        R_alk1 = self.rotationMatrix(templateAlk1Vector, CGThioAlk1)

        for atom in self.monomerAtoms:
            transformedAtom = copy.deepcopy(atom)
            transformedAtomCoords = np.transpose(R_alk1*np.transpose(np.matrix(transformedAtom[1])))
            transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
            alk1ManipulatedMonomer.append(transformedAtom)


        thioNormalVector = self.determineNormalToThiopheneRing(alk1ManipulatedMonomer)
#        thioNormalVector = self.determineNormalToThiopheneRing(alk1ManipulatedMonomer, CGThioAlk1) # We just mapped the templateAlk1Vector to CGThioAlk1, so its Alk1 vector should already be CGThioAlk1

        R_thio = self.rotationMatrix(thioNormalVector, CGThioPlaneNormal)

        # Try only rotating the thiophene ring and leave the alkyl sidechain in its place
        thioRingElements = [0, 1, 3, 4, 5, 7] # These are the line numbers of the atoms that belong to the thiophene ring
        alk1Elements = [8, 2, 12, 10, 6, 15, 13, 9, 17] # These line numbers are the alk1 atoms
        alk2Elements = [16, 11, 20, 19, 14, 22, 21, 18, 23, 24] # These line numbers are the alk2 atoms

        for atomNo in range(len(alk1ManipulatedMonomer)):
            transformedAtom = copy.deepcopy(alk1ManipulatedMonomer[atomNo])
            transformedAtomCoords = np.transpose(R_thio*np.transpose(np.matrix(transformedAtom[1])))
            transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
            alk1AndThioManipulatedMonomer.append(transformedAtom)

        # # ORCA inputs still don't work, so now we need to rotate the ring around the normal to this ring

        newThioNormalVector = self.determineNormalToThiopheneRing(alk1AndThioManipulatedMonomer)
        
        # for atom in alk1AndThioManipulatedMonomer:
        #     alk1AndThioFullyManipulatedMonomer.append(atom)


        if self.rotateThio == True:

            R_thioZ = self.thioInPlaneRotation(alk1AndThioManipulatedMonomer, CGAdjThioCoords, newThioNormalVector)

            for atomNo in range(len(alk1AndThioManipulatedMonomer)):
                transformedAtom = copy.deepcopy(alk1AndThioManipulatedMonomer[atomNo])
                if atomNo in thioRingElements: #  Just Rotate the Thiophene Ring
                    transformedAtomCoords = np.transpose(R_thioZ*np.transpose(np.matrix(transformedAtom[1])))
                    transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
                alk1AndThioZManipulatedMonomer.append(transformedAtom)
        else:
            for atom in alk1AndThioManipulatedMonomer:
                alk1AndThioZManipulatedMonomer.append(atom)



        if self.rotateAlk1 == True:

            # ---=== How to fix the rest of the geometry: ===---
            # 0) Calculate the current COM position of the Alk1 chain (C9, H3, H13, C11, H7, H16, C14, H10, H18)
            currentAlk1COM = self.determineCOM(alk1AndThioZManipulatedMonomer, alk1Only = True) # As if the origin is the thio
            # 0.5) Move all of the Alk1 atoms so that their new COM position matches that of COMAlk1
            alk1Movement = [(relativeAlk1Posn[0] - currentAlk1COM[0]), (relativeAlk1Posn[1] - currentAlk1COM[1]), (relativeAlk1Posn[2] - currentAlk1COM[2])]
            # 1) Make a coordinate switch so that COMAlk1 is now the origin (Move all atoms st. CurrentCoordinates = CurrentCoordinate - relativeAlk1Posn)
            for atomNo in range(len(alk1AndThioZManipulatedMonomer)):
                if atomNo in alk1Elements:
                    alk1AndThioZManipulatedMonomer[atomNo][1] = list(np.array(alk1AndThioZManipulatedMonomer[atomNo][1]) + np.array(alk1Movement))
                alk1AndThioZManipulatedMonomer[atomNo][1] = list(np.array(alk1AndThioZManipulatedMonomer[atomNo][1]) - np.array(relativeAlk1Posn))
            # 2) "CurrentVector" = findAxis( C14, C9 )
            currentAlk1Vector = findAxis(alk1AndThioZManipulatedMonomer[13][1], alk1AndThioZManipulatedMonomer[8][1])
            # 3) "TargetVector" = findAxis( COMAlk1, C6 ) # NB COMAlk1 == 0
            # Might not be C6, sometimes ThioZ flips the pentagon but as it is regular, we can just find the closest carbon out of C5, C6
            C5Sep = calcSeparation(relativeAlk1Posn, alk1AndThioZManipulatedMonomer[4][1])
            C6Sep = calcSeparation(relativeAlk1Posn, alk1AndThioZManipulatedMonomer[5][1])
            if (C5Sep < C6Sep):
                # Connect to C5
                targetVector = findAxis(np.array([0,0,0]), alk1AndThioZManipulatedMonomer[4][1])
                # targetVector = findAxis(np.array(relativeAlk1Posn), alk1AndThioZManipulatedMonomer[4][1])
            else:
                # Connect to C6
                targetVector = findAxis(np.array([0,0,0]), alk1AndThioZManipulatedMonomer[5][1])
                # targetVector = findAxis(np.array(relativeAlk1Posn), alk1AndThioZManipulatedMonomer[5][1])
            # 4) Map currentVector onto targetVector as before
            R_alk1_rot = self.rotationMatrix(currentAlk1Vector, targetVector)

            for atomNo in range(len(alk1AndThioZManipulatedMonomer)):
                transformedAtom = copy.deepcopy(alk1AndThioZManipulatedMonomer[atomNo])
                if atomNo in alk1Elements:
                    transformedAtomCoords = np.transpose(R_alk1_rot*np.transpose(np.matrix(transformedAtom[1])))
                    transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
                # 5) Move all atoms by relativeAlk1Posn (current += rel) to return the ThioCOM to the origin again
                transformedAtom[1] = list(np.array(transformedAtom[1]) + np.array(relativeAlk1Posn))
                thioAlk1CompletedMonomer.append(transformedAtom)

        else:
            for atom in alk1AndThioZManipulatedMonomer:
                thioAlk1CompletedMonomer.append(atom)



        if self.rotateAlk2 == True:


            # A) REPEAT ALL STEPS AGAIN FOR THE ALK2 CHAIN
            # 0) Calculate the current COM position of the Alk2 chain = ( C17, H12, H21, C20, H15, H23, C22, H19, H24, H25 )
            currentAlk2COM = self.determineCOM(thioAlk1CompletedMonomer, alk2Only = True) # As if the origin is the thio
            # 0.5) Move all of the Alk2 atoms so that their new COM position matches that of COMAlk1
            alk2Movement = [(relativeAlk2Posn[0] - currentAlk2COM[0]), (relativeAlk2Posn[1] - currentAlk2COM[1]), (relativeAlk2Posn[2] - currentAlk2COM[2])]
            # 1) Make a coordinate switch so that COMAlk2 is now the origin (Move all atoms st. CurrentCoordinates = CurrentCoordinate - relativeAlk2Posn)
            for atomNo in range(len(thioAlk1CompletedMonomer)):
                if atomNo in alk2Elements:
                    thioAlk1CompletedMonomer[atomNo][1] = list(np.array(thioAlk1CompletedMonomer[atomNo][1]) + np.array(alk2Movement))
                thioAlk1CompletedMonomer[atomNo][1] = list(np.array(thioAlk1CompletedMonomer[atomNo][1]) - np.array(relativeAlk2Posn))
            # 2) "CurrentVector" = findAxis( C22, C17 )
            currentAlk2Vector = findAxis(thioAlk1CompletedMonomer[21][1], thioAlk1CompletedMonomer[16][1])
            # 3) "TargetVector" = findAxis( COMAlk2, C14 ) # NB COMAlk2 == 0
            targetVector = findAxis(np.array([0, 0, 0]), thioAlk1CompletedMonomer[13][1])

            # targetVector = findAxis(np.array(relativeAlk2Posn), thioAlk1CompletedMonomer[13][1])
            # 4) Map currentVector onto targetVector as before
            R_alk2_rot = self.rotationMatrix(currentAlk2Vector, targetVector)

            for atomNo in range(len(thioAlk1CompletedMonomer)):
                transformedAtom = copy.deepcopy(thioAlk1CompletedMonomer[atomNo])
                if atomNo in alk2Elements:
                    transformedAtomCoords = np.transpose(R_alk2_rot*np.transpose(np.matrix(transformedAtom[1])))
                    transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
                # 5) Move all atoms by relativeAlk2Posn (current += rel) to return the ThioCOM to the origin again
                transformedAtom[1] = list(np.array(transformedAtom[1]) + np.array(relativeAlk2Posn))
                thioAlk1Alk2CompletedMonomer.append(transformedAtom)

        else:
            for atom in thioAlk1CompletedMonomer:
                thioAlk1Alk2CompletedMonomer.append(atom)


        # As a final check, somehow after the thioZ rotation, the hydrogen atom ends up being bonded to the wrong carbon and is too close to the alk1 chain.
        # This probably means that the thioZ is rotating the wrong way, but as we're dealing with a regular pentagon, we could quite easily just flip the hydrogen round (as it also makes basically no contribution to the ZINDO measurements, but could possibly prevent the SCF converging if it's located within another atom)
        # To this end, rotate the hydrogen atom (H8) 180 degrees around the axis S1 -> <C6, C5> if the separation between H8 and C9 is < 1A
        separationCheck = calcSeparation(thioAlk1Alk2CompletedMonomer[7][1], thioAlk1Alk2CompletedMonomer[8][1])
        if separationCheck < 1:
            sulphurAtom = thioAlk1Alk2CompletedMonomer[0][1]
            # Opposite side is the average position of C5 and C6
            oppositeSide = [(thioAlk1Alk2CompletedMonomer[4][1][0] + thioAlk1Alk2CompletedMonomer[5][1][0])/2., (thioAlk1Alk2CompletedMonomer[4][1][1] + thioAlk1Alk2CompletedMonomer[5][1][1])/2., (thioAlk1Alk2CompletedMonomer[4][1][2] + thioAlk1Alk2CompletedMonomer[5][1][2])/2.]
            thioSplitAxis = findAxis(sulphurAtom, oppositeSide)
            R_H8_rot = self.rotationMatrixAroundAxis(np.pi, thioSplitAxis)
            transformedAtom = copy.deepcopy(thioAlk1Alk2CompletedMonomer[7])
            transformedAtomCoords = np.transpose(R_H8_rot*np.transpose(np.matrix(transformedAtom[1])))
            transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
            thioAlk1Alk2CompletedMonomer[7] = transformedAtom


        # Add on an additional Hydrogen Atom
        
        if (self.justPlotThios == True):

            tempHydrogen = copy.deepcopy(thioAlk1Alk2CompletedMonomer[7]) #  No longer atom[7] - it's still H8, but we cut out atoms so it's now just the last one
            sulphurAtom = thioAlk1Alk2CompletedMonomer[0][1]
            oppositeSide = [(thioAlk1Alk2CompletedMonomer[4][1][0] + thioAlk1Alk2CompletedMonomer[5][1][0])/2., (thioAlk1Alk2CompletedMonomer[4][1][1] + thioAlk1Alk2CompletedMonomer[5][1][1])/2., (thioAlk1Alk2CompletedMonomer[4][1][2] + thioAlk1Alk2CompletedMonomer[5][1][2])/2.]
            thioSplitAxis = findAxis(sulphurAtom, oppositeSide)
            R_H8_rot = self.rotationMatrixAroundAxis(np.pi, thioSplitAxis)
            tempHydrogenCoords = np.transpose(R_H8_rot*np.transpose(np.matrix(tempHydrogen[1])))
            tempHydrogen[1] = [tempHydrogenCoords[0,0], tempHydrogenCoords[0,1], tempHydrogenCoords[0,2]]
            thioAlk1Alk2CompletedMonomer.append(tempHydrogen)




        if (self.justPlotThios == True):
            for atomNo in range(len(thioAlk1Alk2CompletedMonomer)):
                if atomNo in thioRingElements: # Just the thiophene ring
                    monomerFinal.append(thioAlk1Alk2CompletedMonomer[atomNo])
            monomerFinal.append(thioAlk1Alk2CompletedMonomer[-1]) # The extra Hydrogen
        else:
            for atom in thioAlk1Alk2CompletedMonomer:
                monomerFinal.append(atom)
        # for atomNo in range(len(thioAlk1Alk2CompletedMonomer)):
        #     # if (atomNo == 7):
        #     #     continue
        #     # if (atomNo in thioRingElements):
        #     testMonomer.append(thioAlk1Alk2CompletedMonomer[atomNo])


        # return testMonomer




        # Finally add on the ThioCOM Coords

        for atom in monomerFinal:
            atom[1] = [atom[1][0]+COMThio[0], atom[1][1]+COMThio[1], atom[1][2]+COMThio[2]]





        return monomerFinal



# --------===============================-------





    # def manipulateMonomer(self, COMCoord, CGNormal, CGThioAlk1):
    #     torsionManipulatedMonomer = [] # The first rotation stage
    #     alk1ManipulatedMonomer = [] # The second rotation stage

 
    #     # Sort out rotation first
    #     templateNormal = self.determineNormalToThiopheneRing(self.monomerAtoms)

    #     templateAlk1Vector = self.determineThioToAlk1Vector(self.monomerAtoms)


    #     print "Template angle between normal and Alk1 =", np.dot(templateNormal, templateAlk1Vector)
    #     print "CG angle between normal and Alk1 =", np.dot(CGNormal, CGThioAlk1)



    #     raise SystemError('STOP')

    #     # Sorting out rotation is none trivial - used this: http://goo.gl/7kIAEp
    #     # v = a x b, therefore s = sin(theta) = |v| and c = cos(theta) = a.b
    #     # Rotation matrix R = I + [V] + [V]**2 * (1 - c)/s**2
    #     # Where [V] is the "skew-symmetric cross-product matrix of v" =
    #     # [ [ 0, -v_{3}, v_{2} ], [ v_{3}, 0, -v_{1} ], [ -v_{2}, v_{1}, 0 ] ]

    #     # print "TemplateNormal =", templateNormal
    #     # print "CGNormal =", CGNormal

    #     crossProduct = np.cross(templateNormal, CGNormal)
    #     sinAngle = np.sqrt(((crossProduct[0])**2) + ((crossProduct[1])**2) + ((crossProduct[2])**2))
    #     cosAngle = np.dot(templateNormal, CGNormal)

    #     skewMatrix = np.matrix([[0, -crossProduct[2], crossProduct[1]], [crossProduct[2], 0, -crossProduct[0]], [-crossProduct[1], crossProduct[0], 0]])
    #     skewMatrixSquared = skewMatrix*skewMatrix

    #     R_Torsion = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix + skewMatrixSquared*((1 - cosAngle)/(sinAngle**2))

    #     # print "cosAngle =", cosAngle
    #     # print "sinAngle =", sinAngle
    #     # print "prefactor =", (1 - cosAngle)/(sinAngle**2)
    #     # print "Skew =", skewMatrix
    #     # print "Skew Squared =", skewMatrixSquared
    #     # print "Prefactor * skeqSquared =", skewMatrixSquared*((1-cosAngle)/(sinAngle**2))
    #     # print "1+d =", np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix


    #     # print "---======---"

    #     # print "R =", R_Torsion

    #     # raise SystemError('STOP')

    #     # Apply the torsional angle rotation first

    #     for atom in self.monomerAtoms:
    #         transformedAtom = copy.deepcopy(atom)
    #         transformedAtomCoords = (np.transpose(R_Torsion*np.transpose(np.matrix(transformedAtom[1]))))
    #         # print "inital Coords =", transformedAtom[1]
    #         # print "transformed Coords =", transformedAtomCoords
    #         transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
    #         # print "COM =", COMCoord
    #         # print "Final =", transformedAtom[1]
    #         # raise SystemError('STOP')
    #         torsionManipulatedMonomer.append(transformedAtom)
    #         # print "Atom =", atom, "Manipulated =", transformedAtom


    #     # NOW need to find out where the thio to Alk vector is pointing for the manipulated monomer

    #     templateAlk1Vector = self.determineThioToAlk1Vector(torsionManipulatedMonomer)

    #     # Now calculate the rotation matrix to rotate the alkyl vector into place

    #     # print "templateAlk1Vector =", templateAlk1Vector
    #     # print "CGThioAlk1 =", CGThioAlk1

    #     # print "torsionManipulatedMonomer =", torsionManipulatedMonomer
        
    #     crossProduct = np.cross(templateAlk1Vector, CGThioAlk1)
    #     sinAngle =  np.sqrt(((crossProduct[0])**2) + ((crossProduct[1])**2) + ((crossProduct[2])**2))
    #     cosAngle = np.dot(templateAlk1Vector, CGThioAlk1)

    #     skewMatrix = np.matrix([[0, -crossProduct[2], crossProduct[1]], [crossProduct[2], 0, -crossProduct[0]], [-crossProduct[1], crossProduct[0], 0]])
    #     skewMatrixSquared = skewMatrix*skewMatrix

    #     R_Alk1 = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix + skewMatrixSquared*((1 - cosAngle)/(sinAngle**2))

    #     # Now apply the alk1 rotation into place and add on the COM of the CG Coords

    #     for atom in torsionManipulatedMonomer:
    #         transformedAtom = copy.deepcopy(atom)
    #         transformedAtomCoords = (np.transpose(R_Alk1*np.transpose(np.matrix(transformedAtom[1]))))
    #         transformedAtom[1] = [transformedAtomCoords[0,0]+COMCoord[0], transformedAtomCoords[0,1]+COMCoord[1], transformedAtomCoords[0,2]+COMCoord[2]]
    #         alk1ManipulatedMonomer.append(transformedAtom)


    #     # print "alk1ManipulatedMonomer =", alk1ManipulatedMonomer


    #     # print "Target Thio Normal =", CGNormal
    #     # print "Actual Thio Normal =", self.determineNormalToThiopheneRing(alk1ManipulatedMonomer)

    #     # print "Target Alk1 Vector =", CGThioAlk1
    #     # print "Actual Alk1 Vector =", self.determineThioToAlk1Vector(alk1ManipulatedMonomer)





    #     # # print "Self.MonomerAtoms =", self.monomerAtoms
    #     # # print "ManipulatedMonomer =", manipulatedMonomer

    #     # # currentAlk1Vector = self.determineThioToAlk1Vector(manipulatedMonomer)

    #     # # print "Current Alk1 Vector =", currentAlk1Vector, np.sqrt((currentAlk1Vector[0]**2)+(currentAlk1Vector[1]**2)+(currentAlk1Vector[2]**2))
    #     # # print "Target Alk1 Vector =", CGThioAlk1, np.sqrt((CGThioAlk1[0]**2)+(CGThioAlk1[1]**2)+(CGThioAlk1[2]**2))
    #     # # print "Normal vector axis =", CGNormal, np.sqrt((CGNormal[0]**2)+(CGNormal[1]**2)+(CGNormal[2]**2))

    #     # # raise SystemError('STO}')




    #     # # and add on the CGCoord to each atom in self.monomerAtoms

    #     return alk1ManipulatedMonomer
        


    def generateORCAInput(self):
        # CreateName (rounded to .1 Ang)
        self.createName()
        # Check that file with the right name doesn't already exist
        #       If it does, pass.
        #       Otherwise make the ORCA inputFile
        exists = self.makeDirs()
        if exists == True:
            print "File", self.fullPath, "already exists, skipping...\n"
            return
        atomsToWrite = self.getAtomsToWrite()
        self.writeInpFile(atomsToWrite)
        # RUN ORCA
        # Analyse ORCA outputs to create a structure (dictionary?) that takes two segment numbers and returns the transferIntegral        



    def writeInpFile(self, atomsToWrite):
        templateFile = open('./templates/template.inp', 'r')
        templateLines = templateFile.readlines()
        templateFile.close()
        linesToWrite = []
        for lineNo in range(len(templateLines)):
            if templateLines[lineNo] == '*':
                for atom in atomsToWrite:
                    lineToWrite = ' '
                    lineToWrite += str(atom[0])+' '
                    lineToWrite += str(atom[1])+' '
                    lineToWrite += str(atom[2])+' '
                    lineToWrite += str(atom[3])+'\n'
                    linesToWrite.append(lineToWrite)
            linesToWrite.append(templateLines[lineNo])
        orcaInputFile = open(self.fullPath, 'w+')
        orcaInputFile.writelines(linesToWrite)
        orcaInputFile.close()
        print "Orca Input file written to", self.fullPath, "\n"
        


    def getAtomsToWrite(self):
        atomsToWrite = []
        segmentNo = 1
        atomNo = 1
        for segment in self.segmentAtoms:
            for CGAtom in segment.iteritems():
                # print "Treating segment", str(segmentNo)+", monomer", str(atomNo)+"."
                CGCoord = np.array([CGAtom[1][0], CGAtom[1][1], CGAtom[1][2]])
                CGThioPlaneNormal = CGAtom[1][5]
                CGThioAlk1Axis = CGAtom[1][6]
                CGAdjThioCoords = CGAtom[1][7]
                CGAlk1 = CGAtom[1][8]
                CGAlk2 = CGAtom[1][9]
                fineGrainedAtoms = self.manipulateMonomer(CGCoord, CGThioPlaneNormal, CGAdjThioCoords, CGThioAlk1Axis, CGAlk1, CGAlk2)
                for atom in fineGrainedAtoms:
                    atomsToWrite.append([atom[0], str(atom[1][0]), str(atom[1][1]), str(atom[1][2])])
                atomNo += 1
            segmentNo += 1
        # Write ORCA File
        return atomsToWrite



def plotTest(segmentsMaster, periodicSegmentsMaster):
    segment5 = []
    segment266 = []
    segment4 = []
    segment112 = []
    segment113 = []
    segment104 = []
    periodicSegment112 = []
    periodicSegment113 = []
    periodicSegment104 = []
    for segment in segmentsMaster:
        for atom in segment.iteritems():
            if atom[1][4] == 5:
                segment5.append([atom[1][0], atom[1][1], atom[1][2]])
            elif atom[1][4] == 266:
                segment266.append([atom[1][0], atom[1][1], atom[1][2]])
            elif atom[1][4] == 4:
                segment4.append([atom[1][0], atom[1][1], atom[1][2]])
            elif atom[1][4] == 112:
                segment112.append([atom[1][0], atom[1][1], atom[1][2]])
            elif atom[1][4] == 113:
                segment113.append([atom[1][0], atom[1][1], atom[1][2]])
            elif atom[1][4] == 104:
                segment104.append([atom[1][0], atom[1][1], atom[1][2]])
    for segment in periodicSegmentsMaster.items():
        if (segment[0] == 112) or (segment[0] == 113) or (segment[0] == 104):
            for atom in segment[1].items():
                if atom[1][4] == 112:
                    periodicSegment112.append([atom[1][0], atom[1][1], atom[1][2]])
                elif atom[1][4] == 113:
                    periodicSegment113.append([atom[1][0], atom[1][1], atom[1][2]])
                elif atom[1][4] == 104:
                    periodicSegment104.append([atom[1][0], atom[1][1], atom[1][2]])
    # SEGMENT 5
    redX = []
    redY = []
    redZ = []
    # NEARBY SEGMENTS
    blueX = []
    blueY = []
    blueZ = []
    # REAL COUNTERPART OF PERIODICS
    greenX = []
    greenY = []
    greenZ = []
    # PERIODICS
    blackX = []
    blackY = []
    blackZ = []
    

    for atom in segment5:
        redX.append(atom[0])
        redY.append(atom[1])
        redZ.append(atom[2])
        
    for atom in segment4+segment266:
        blueX.append(atom[0])
        blueY.append(atom[1])
        blueZ.append(atom[2])

    for atom in segment112+segment113+segment104:
        greenX.append(atom[0])
        greenY.append(atom[1])
        greenZ.append(atom[2])

    for atom in periodicSegment112+periodicSegment113+periodicSegment104:
        blackX.append(atom[0])
        blackY.append(atom[1])
        blackZ.append(atom[2])
        

    fig = P.figure()
    ax = p3.Axes3D(fig)
    ax.scatter(redX, redY, redZ, s = 80, c = 'r')
    ax.scatter(greenX, greenY, greenZ, s = 40, c = 'g')
    ax.scatter(blueX, blueY, blueZ, s = 40, c = 'b')
    ax.scatter(blackX, blackY, blackZ, s = 40, c = 'k')

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    P.savefig('./test.png')

    




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dat", help="The data file *.dat to run the paracrystallinity measurements on")
    parser.add_argument("-l", "--lammpstrj", help="The corresponding *.lammpstrj file that includes all of the atomic trajectories to run the paracrystallinity measurements on")
    args = parser.parse_args()

    testingIndividualMolecules = False
    testMoleculeNumber = 1

    rotateThio = True
    rotateAlk1 = False
    rotateAlk2 = False
    justPlotThios = True

    # toleranceAngles = np.arange(np.pi/180., 46*np.pi/180., np.pi/180)

    toleranceAngle = np.pi/6.
    maximumHoppingDistance = 15


    ###########################################################################################
    #### THIS SECTION OF CODE RAN ALL OF THE LAMMPSTRJ FILES IN THIS DIRECTORY IN SEQUENCE ####
    ###########################################################################################
    # datFiles, lammpstrjFiles = getFilesList('./')
    # print "Valid files =", zip(datFiles, lammpstrjFiles)

    # for morphologyNo in range(len(datFiles)):

    #     datName = datFiles[morphologyNo]
    #     lammpstrjName = lammpstrjFiles[morphologyNo]
    ###########################################################################################


    datName = args.dat
    lammpstrjName = args.lammpstrj

    segmentsMaster = [] # This should just have all of the segments in the morphology, vaguely sorted by molecule but not split by molecule
    atomMaster = {}

    print "Running segments.py for", datName, "and", lammpstrjName

    if testingIndividualMolecules == True:
        trimTrajectory(lammpstrjName, testMoleculeNumber)
        trajectoryData, allTimestepNos, simVolData = loadTrajectory('./singleMolTraj.lammpstrj')

    else:
        trajectoryData, allTimestepNos, simVolData = loadTrajectory(lammpstrjName)

    masses, atoms, bonds, angles, dihedrals, impropers = loadDat(datName)
    thiosOnly = getAtomTypes(trajectoryData[-1], 'thio') # Take final timestep
    alk1sOnly = getAtomTypes(trajectoryData[-1], 'alk1')
    alk2sOnly = getAtomTypes(trajectoryData[-1], 'alk2')
    # Results in thiophenes or alk1s sorted by molecule

    globalAverageSegmentLength = []
    globalCoordinatePositions = {}
    rollingSegmentNumber = 1

    for moleculeThiosNo in range(len(thiosOnly)):
        print "Treating molecule", str(moleculeThiosNo+1) + "..."

        thisMolecule = molecule(thiosOnly[moleculeThiosNo], alk1sOnly[moleculeThiosNo], alk2sOnly[moleculeThiosNo], bonds, toleranceAngle)

        segments = thisMolecule.returnSegments()
        thioOrientations, alkOrientations, neighbouringThioCoords, bondedAlk1Coords, bondedAlk2Coords = thisMolecule.returnOrientations()

        segmentLengths = []
        for segment in segments:
            for atom in segment.iteritems():
                segment[atom[0]].append(moleculeThiosNo+1)
                segment[atom[0]].append(rollingSegmentNumber)
                segment[atom[0]].append(thioOrientations[atom[0]])
                segment[atom[0]].append(alkOrientations[atom[0]])
                segment[atom[0]].append(neighbouringThioCoords[atom[0]])
                segment[atom[0]].append(bondedAlk1Coords[atom[0]])
                segment[atom[0]].append(bondedAlk2Coords[atom[0]])
                atomMaster[atom[0]] = copy.deepcopy(atom[1])
                globalCoordinatePositions[(atom[1][0], atom[1][1], atom[1][2])] = rollingSegmentNumber
            rollingSegmentNumber += 1
            segmentLengths.append(len(segment))
            segmentsMaster.append(segment)
            csvFile = open('./orcaInputs/'+datName[:-4]+'/SegmentLengths.csv', 'a+')
            csvWriter = csv.writer(csvFile, delimiter = ',')
            csvWriter.writerow([rollingSegmentNumber-1, len(segment)])
            csvFile.close()

    ###### SEGMENTS MASTER FORMAT #####
    # {AtomIndex: [XCoord, YCoord, ZCoord, moleculeNumber, segmentNumber, [VectorBetweenBothAdjacentThios(PlaneOfThiopheneRing)], [VectorFromThioToBondedAlk1], [CoordinatesOfFirstNeighbouringThiophene, CoordinatesOfSecondNeighbouringThiopehen(ifPresent)], [Coordinates of bonded Alk1], [Coordinates of bonded Alk2]]}
    ###################################



        averageSegmentLength = np.average(segmentLengths)
        print "Average segment length for this molecule =", averageSegmentLength

        globalAverageSegmentLength.append(averageSegmentLength)

        fileName = getFilename(args.dat, moleculeThiosNo)

        # plotMolecule(thiosOnly, segments, fileName)

    print "Writing segment coordinates to csv..."
    csvFile = open('./'+datName[:-4]+'_SegmentCoordinates.csv', 'w+')
    csvWriter = csv.writer(csvFile, delimiter = ',')
    csvWriter.writerows(globalCoordinatePositions.items())
    csvFile.close()
    print "Segment coordinates written to ./"+str(datName[:-4])+"_SegmentCoordinates.csv."



    csvFileName = fileName[:-7]+'_tolerances'

    print "\n", globalAverageSegmentLength
    averageSegmentLengthGlobal = np.average(globalAverageSegmentLength)
    print "Average global segment length for this morphology =", averageSegmentLengthGlobal

    # print "Tol angle =", toleranceAngle
    # writeCSV(toleranceAngle, averageSegmentLengthGlobal, csvFileName)

    # ---=== TO DO ===---

    # Need to find pairs of these segments to be submitted into ORCA


    print "There are", len(segmentsMaster), "segments."

    # SegmentsMaster is now a list of dictionaries. Each dictionary in the list represents a new segment. Within each segment, the key is the LAMMPS Atomindex, and the value is the following list: [XCoord, YCoord, ZCoord, MoleculeNo, SegmentNo]
    # csvFileName = fileName[:-7]+'_segmentPairs'
    # totalNumberOfSegmentPairs = []
    # for maximumHoppingDistance in maximumHoppingDistances:
        # totalNumberOfSegmentPairs.append(numberOfSegmentPairs)
    #     writeCSV(maximumHoppingDistance, numberOfSegmentPairs, csvFileName)
    # print totalNumberOfSegmentPairs

    # REMOVING THIS FOR THE TIME BEING AS IT'S SLOW
    t1 = T.time()
    segmentPairs, periodicSegmentPairs, periodicSegmentsMaster = findSegmentPairs(segmentsMaster, atomMaster, maximumHoppingDistance, simVolData[0])
    t2 = T.time()
    print "Took %.1f seconds to find segments." % (t2-t1)

    # segmentPairs = [[1], [1,2], [1,64], [1,333], [1,357], [1,356]]
    # segmentPairs = [[928]]

    # SegementPairs also includes each individual segment too, to make sure that we have all of the data.

    print "\nNow treating segment pairs..."

    t0 = T.time()
    totalPairs = len(segmentPairs)+len(periodicSegmentPairs)
    completedPairs = 0

    for segmentPair in segmentPairs:
        print "Writing input for segment(s):", segmentPair
        t3 = T.time()
        ORCAInput(segmentsMaster, segmentPair, datName, rotateThio, rotateAlk1, rotateAlk2, justPlotThios).generateORCAInput()
        t4 = T.time()
        completedPairs += 1
        currentElapsedTime = float(t4-t0)
        if currentElapsedTime < 60.0:
            timeunits = 'seconds.'
        elif currentElapsedTime < 3600:
            currentElapsedTime /= 60.0
            timeunits = 'minutes.'
        elif currentElapsedTime < 86400:
            currentElapsedTime /= 3600.0
            timeunits = 'hours.'
        else:
            currentElapsedTime /= 86400.0
            timeunits = 'days.'
        averageTimePerPair = float(t4-t0)/float(completedPairs)
        expectedTotalTime = float(averageTimePerPair)*float(totalPairs)
        remainingTime = expectedTotalTime - float(t4-t0)
        if remainingTime < 60.0:
            timeunits2 = 'seconds.'
        elif remainingTime < 3600:
            remainingTime /= 60.0
            timeunits2 = 'minutes.'
        elif remainingTime < 86400:
            remainingTime /= 3600.0
            timeunits2 = 'hours.'
        else:
            remainingTime /= 86400.0
            timeunits2 = 'days.'

        percentageComplete = int((completedPairs/float(totalPairs))*100.0)
        previousOutputtedPercentage = 0
        previousOutputtedPercentage = percentageComplete
        print "---======---"
        print "Calculations %d percent complete. Current elapsed time = %.1f %s Estimated time remaining = %.1f %s" % (percentageComplete, currentElapsedTime, timeunits, remainingTime, timeunits2)
        print "---======---\n"


    print "-----==========-----"
    print "THE REAL SEGMENT PAIRS HAVE NOW BEEN TREATED, NOW IT IS TIME TO TREAT THE PERIODIC SEGMENTS"
    print "-----==========-----"



    for segmentPair in periodicSegmentPairs:
        print "Writing input for segment(s):", segmentPair
        t3 = T.time()
        ORCAInput(segmentsMaster, segmentPair, datName, rotateThio, rotateAlk1, rotateAlk2, justPlotThios, periodicSegmentPair = True, periodicMasterDict = periodicSegmentsMaster).generateORCAInput()
        t4 = T.time()
        completedPairs += 1
        currentElapsedTime = float(t4-t0)
        if currentElapsedTime < 60.0:
            timeunits = 'seconds.'
        elif currentElapsedTime < 3600:
            currentElapsedTime /= 60.0
            timeunits = 'minutes.'
        elif currentElapsedTime < 86400:
            currentElapsedTime /= 3600.0
            timeunits = 'hours.'
        else:
            currentElapsedTime /= 86400.0
            timeunits = 'days.'
        averageTimePerPair = float(t4-t0)/float(completedPairs)
        expectedTotalTime = float(averageTimePerPair)*float(totalPairs)
        remainingTime = expectedTotalTime - float(t4-t0)
        if remainingTime < 60.0:
            timeunits2 = 'seconds.'
        elif remainingTime < 3600:
            remainingTime /= 60.0
            timeunits2 = 'minutes.'
        elif remainingTime < 86400:
            remainingTime /= 3600.0
            timeunits2 = 'hours.'
        else:
            remainingTime /= 86400.0
            timeunits2 = 'days.'

        percentageComplete = int((completedPairs/float(totalPairs))*100.0)
        previousOutputtedPercentage = 0
        previousOutputtedPercentage = percentageComplete
        print "---======---"
        print "Calculations %d percent complete. Current elapsed time = %.1f %s Estimated time remaining = %.1f %s" % (percentageComplete, currentElapsedTime, timeunits, remainingTime, timeunits2)
        print "---======---\n"









    # Need to output these segments in a form that the CT program can understand
