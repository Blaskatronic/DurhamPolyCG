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
    def __init__(self, thiosInThisMolecule, alk1sInThisMolecule,  bonds, toleranceAngle):
        self.allThios = self.makeDict(thiosInThisMolecule)
#        self.thiosInMorphology = self.makeDict(thiosInMorphology)
        self.thios = copy.deepcopy(self.allThios)
        self.alk1s = self.makeDict(alk1sInThisMolecule)
        self.thioThioBonds, self.thioAlk1Bonds = self.findRelevantBonds(bonds)
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
        for bond in bonds:
            if (bond[1] == 1): # Thio-thio bonds only
                if (bond[2] in self.thios) or (bond[3] in self.thios):
                    thioThioBonds.append([bond[2], bond[3]])
            elif (bond[1] == 2): # Thio-Alk1 bonds only
                if (bond[2] in self.thios):
                    thioAlkBonds.append([bond[2], bond[3]])
                elif (bond[3] in self.thios):
                    thioAlkBonds.append([bond[3], bond[2]])
        return thioThioBonds, thioAlkBonds

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
        thioPlaneOrientationMaster = {}
        alkylSideChainOrientationMaster = {}
        neighbouringThiosCoords = {}
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
                        break
                    elif atomNo == bond[1]:
                        alk1Coords = self.alk1s[bond[0]]
                        break

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
        return thioPlaneOrientationMaster, alkylSideChainOrientationMaster, neighbouringThiosCoords



    def returnSegments(self):
        segmentReturn = []
        for segment in self.segments:
            segmentReturn.append({})
            for atom in segment:
                segmentReturn[-1][atom] = self.allThios[atom]
        return segmentReturn

    def returnOrientations(self):
        thioOrientations, alkylSideChainOrientations, neighbouringThiosCoords = self.findOrientations()
        return thioOrientations, alkylSideChainOrientations, neighbouringThiosCoords


        
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



def findSegmentPairs(segmentsMaster, atomMaster, maximumHoppingDistance):
    segmentPairs = []
    # Just look at first one first
    for i in range(len(segmentsMaster)):
        segment = segmentsMaster[i]
        nearbySegments = []

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

        print "Current Segment =", currentSegment, "Nearby Segments =", nearbySegments

    print "SegmentPairs =", segmentPairs
    print "Len SegmentPairs =", len(segmentPairs)
    print "Job done."
    return segmentPairs



# -------------==================== ======================-----------------




# -------------==================== Creating ORCA Inputs  ======================-----------------

class ORCAInput:
    def __init__(self, segmentsMaster, segmentPair):
        self.segmentPair = segmentPair
        self.segmentAtoms = []
        self.segmentAtoms.append(segmentsMaster[segmentPair[0]-1])
        self.segmentAtoms.append(segmentsMaster[segmentPair[1]-1])
        self.monomerAtoms = self.read3HTTemplate() # These are the coordinates of the atoms representing a CG site at [0, 0, 0]
        self.generateORCAInput()
        

    def returnAnswer(self):
        raise SystemError('DONE FIRST SEGMENT PAIR')
        return 0
        
    def determineCOM(self, inputAtoms, thioOnly = False, alk1Only = False):
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
        # Want a name kinda like C1_LengthOfChain1_C2_LengthOfChain2_X_XSeparationFromCentrePoint_Y_YSep_Z_ZSep.inp
        C1Length = len(self.segmentAtoms[0])
        C2Length = len(self.segmentAtoms[1])
        # Determine COMs of each chain
        chain1 = []
        chain2 = []
        for atom in self.segmentAtoms[0].iteritems():
            chain1.append([0, [atom[1][0], atom[1][1], atom[1][2]]])
        for atom in self.segmentAtoms[1].iteritems():
            chain2.append([0, [atom[1][0], atom[1][1], atom[1][2]]])
        COM1 = self.determineCOM(chain1)
        COM2 = self.determineCOM(chain2)
        separation = COM2-COM1
        self.name = 'C1_'+str(len(chain1))+'_C2_'+str(len(chain2))+'_X_%.1f_Y_%.1f_Z_%.1f' % (float(separation[0]), float(separation[1]), float(separation[2]))


    def checkFileExists(self, name):
        files = os.listdir('./lookupTable')
        if name in files:
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

    def getRotationMatrix(self, theta, axis):
        # Rotation matrix obtained from http://goo.gl/RkW80
        R = np.matrix([[ np.cos(theta)+(axis[0]**2)*(1-np.cos(theta)), axis[0]*axis[1]*(1-np.cos(theta))-axis[2]*np.sin(theta), axis[0]*axis[2]*(1-np.cos(theta))+axis[1]*np.sin(theta) ],
                       [ axis[1]*axis[0]*(1-np.cos(theta))+axis[2]*np.sin(theta), np.cos(theta)+(axis[1]**2)*(1-np.cos(theta)), axis[1]*axis[2]*(1-np.cos(theta))-axis[0]*np.sin(theta) ],
                       [ axis[2]*axis[0]*(1-np.cos(theta))-axis[1]*np.sin(theta), axis[2]*axis[1]*(1-np.cos(theta))+axis[0]*np.sin(theta), np.cos(theta)+(axis[2]**2)*(1-np.cos(theta)) ]])
        return R


    def findRotationAngle(self, atom1, atom2, axis):
        # No idea if this will work, but I got it from http://goo.gl/EVUWsb
        atom1 = np.array(atom1)
        atom2 = np.array(atom2)
        axis = np.array(axis)

        uy = normaliseVec(np.cross(atom1, axis))
        ux = np.cross(uy, axis)

        # This gives us a coordinate triplet system of ux, uy and n which is our new coordinate basis.

        # Atom 1s coordinates are now (ax, 0) = (np.dot(atom1, ux), 0) in the new basis
        # Need to find Atom 2s coordinates:

        bx = np.dot(atom2, ux)
        by = np.dot(atom2, uy)

        theta = math.atan2(by, bx)

        return theta
        

    def thioInPlaneRotation(self, atoms, adjThioCoords, rotationAxis):
#        testRotationAngle = np.pi/180.
        # There are two carbon atoms that connect to adjacent thiophenes:
        leftCarbon = copy.deepcopy(atoms[3])
        rightCarbon = copy.deepcopy(atoms[1])
        # Left and right, however is aribtrary. Therefore, find out which carbon is closest to which adjThio
        if (len(adjThioCoords) == 1) or (len(adjThioCoords) == 2): # For now, just match up one of them.
            # End of chain, therefore line one of them up
            leftSep = calcSeparation(leftCarbon[1], adjThioCoords[0])
            rightSep = calcSeparation(rightCarbon[1], adjThioCoords[0])
            if leftSep <= rightSep:
                # Left carbon is closer to the adjacent thiophene, so find the rotation matrix to minimise this distance
                theta = self.findRotationAngle(leftCarbon[1], adjThioCoords[0], rotationAxis)
                rotationMatrix = self.getRotationMatrix(theta, rotationAxis)
            else:
                theta = self.findRotationAngle(rightCarbon[1], adjThioCoords[0], rotationAxis)
                rotationMatrix = self.getRotationMatrix(theta, rotationAxis)

        # elif len(adjThioCoords) == 2:
        #     pass

        else:
            print "AdjThioCoords =", adjThioCoords
            raise SystemError('Length of AdjThioCoords > 2 or == 0!')

        return rotationMatrix

    def rotationMatrix(self, vector1, vector2):
        # A function to return the rotation matrix around the origin that maps vector1 to vector 2
        crossProduct = np.cross(vector1, vector2)
        sinAngle = np.sqrt(((crossProduct[0]**2) + ((crossProduct[1])**2) + ((crossProduct[2])**2)))
        cosAngle = np.dot(vector1, vector2)

        skewMatrix = np.matrix([[0, -crossProduct[2], crossProduct[1]], [crossProduct[2], 0, -crossProduct[0]], [-crossProduct[1], crossProduct[0], 0]])
        skewMatrixSquared = skewMatrix * skewMatrix
        
        rotMatrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix + skewMatrixSquared*((1 - cosAngle)/(sinAngle**2))

        return rotMatrix



    def manipulateMonomer(self, COMCoord, CGThioPlaneNormal, CGAdjThioCoords, CGThioAlk1):
        alk1ManipulatedMonomer = []
        alk1AndThioManipulatedMonomer = []
        alk1AndThioFullyManipulatedMonomer = []

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

        for atomNo in range(len(alk1ManipulatedMonomer)):
            transformedAtom = copy.deepcopy(alk1ManipulatedMonomer[atomNo])
            transformedAtomCoords = np.transpose(R_thio*np.transpose(np.matrix(transformedAtom[1])))
            transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
            alk1AndThioManipulatedMonomer.append(transformedAtom)

        # # ORCA inputs still don't work, so now we need to rotate the ring around the normal to this ring
        # R_thioZ = self.thioInPlaneRotation(alk1AndThioManipulatedMonomer, CGAdjThioCoords, CGThioPlaneNormal)

        # for atomNo in range(len(alk1ManipulatedMonomer)):
        #     transformedAtom = copy.deepcopy(alk1ManipulatedMonomer[atomNo])
        #     if atomNo in thioRingElements: # Just rotate the Thiophene ring
        #         transformedAtomCoords = np.transpose(R_thioZ*np.transpose(np.matrix(transformedAtom[1])))
        #         transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
        #     alk1AndThioFullyManipulatedMonomer.append(transformedAtom)


        newThioNormalVector = self.determineNormalToThiopheneRing(alk1AndThioManipulatedMonomer)
        
        # for atom in alk1AndThioManipulatedMonomer:
        #     alk1AndThioFullyManipulatedMonomer.append(atom)


        


        # Now add on the COM Coords

        for atom in alk1AndThioFullyManipulatedMonomer:
            atom[1] = [atom[1][0]+COMCoord[0], atom[1][1]+COMCoord[1], atom[1][2]+COMCoord[2]]

        # # I suspect the rotation isn't working properly, therefore insert a few oxygen atoms to show where the normal to the thiophene planes is
        # for i in np.arange(0,3.1,0.2):
        #     alk1AndThioFullyManipulatedMonomer.append(['O', list(np.array(COMCoord)+i*np.array(CGThioPlaneNormal))])


        return alk1AndThioFullyManipulatedMonomer





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
        exists = self.checkFileExists(self.name)
        if exists == True:
            return
        atomsToWrite = self.getAtomsToWrite()
        self.writeInpFile(atomsToWrite)
        # RUN ORCA
        # Analyse ORCA outputs to create a structure (dictionary?) that takes two segment numbers and returns the transferIntegral        
        pass


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
        orcaInputFile = open('./'+str(self.name)+'.inp', 'w+')
        orcaInputFile.writelines(linesToWrite)
        orcaInputFile.close()
        
        


    def getAtomsToWrite(self):
        atomsToWrite = []
        for segment in self.segmentAtoms:
            for CGAtom in segment.iteritems():
                CGCoord = np.array([CGAtom[1][0], CGAtom[1][1], CGAtom[1][2]])
                CGThioPlaneNormal = CGAtom[1][5]
                CGThioAlk1Axis = CGAtom[1][6]
                CGAdjThioCoords = CGAtom[1][7]
                fineGrainedAtoms = self.manipulateMonomer(CGCoord, CGThioPlaneNormal, CGAdjThioCoords, CGThioAlk1Axis)
                for atom in fineGrainedAtoms:
                    atomsToWrite.append([atom[0], str(atom[1][0]), str(atom[1][1]), str(atom[1][2])])
        # Write ORCA File
        return atomsToWrite





if __name__ == '__main__':

    
    testingIndividualMolecules = False
    testMoleculeNumber = 1

    # toleranceAngles = np.arange(np.pi/180., 46*np.pi/180., np.pi/180)

    toleranceAngle = np.pi/7.
    maximumHoppingDistance = 7

    datFiles, lammpstrjFiles = getFilesList('./')
    print "Valid files =", zip(datFiles, lammpstrjFiles)

    for trajectoryNo in range(len(datFiles)):

        datName = datFiles[trajectoryNo]
        lammpstrjName = lammpstrjFiles[trajectoryNo]

        segmentsMaster = [] # This should just have all of the segments in the morphology, vaguely sorted by molecule but not split by molecule
        atomMaster = {}

        print "Treating", datName, "and", lammpstrjName

        if testingIndividualMolecules == True:
            trimTrajectory(lammpstrjName, testMoleculeNumber)
            trajectoryData, allTimestepNos, simVolData = loadTrajectory('./singleMolTraj.lammpstrj')

        else:
            trajectoryData, allTimestepNos, simVolData = loadTrajectory(lammpstrjName)
        
        masses, atoms, bonds, angles, dihedrals, impropers = loadDat(datName)
        thiosOnly = getAtomTypes(trajectoryData[-1], 'thio') # Take final timestep
        alk1sOnly = getAtomTypes(trajectoryData[-1], 'alk1')
        # Results in thiophenes or alk1s sorted by molecule

        globalAverageSegmentLength = []
        rollingSegmentNumber = 1

        for moleculeThiosNo in range(len(thiosOnly)):
            print "Treating molecule", str(moleculeThiosNo+1) + "..."

            thisMolecule = molecule(thiosOnly[moleculeThiosNo], alk1sOnly[moleculeThiosNo], bonds, toleranceAngle)

            segments = thisMolecule.returnSegments()
            thioOrientations, alkOrientations, neighbouringThioCoords = thisMolecule.returnOrientations()

            segmentLengths = []
            for segment in segments:
                for atom in segment.iteritems():
                    segment[atom[0]].append(moleculeThiosNo+1)
                    segment[atom[0]].append(rollingSegmentNumber)
                    segment[atom[0]].append(thioOrientations[atom[0]])
                    segment[atom[0]].append(alkOrientations[atom[0]])
                    segment[atom[0]].append(neighbouringThioCoords[atom[0]])
                    atomMaster[atom[0]] = copy.deepcopy(atom[1])
                rollingSegmentNumber += 1
                segmentLengths.append(len(segment))
                segmentsMaster.append(segment)

                
        ###### SEGMENTS MASTER FORMAT #####
        # {AtomIndex: [XCoord, YCoord, ZCoord, moleculeNumber, segmentNumber, [[VectorBetweenBothAdjacentThios(PlaneOfThiopheneRing)], [VectorFromThioToBondedAlk1], [CoordinatesOfAnyNeighbouringThiophenes(1or2)]]]}
        ###################################



            averageSegmentLength = np.average(segmentLengths)
            print "Average segment length for this molecule =", averageSegmentLength

            globalAverageSegmentLength.append(averageSegmentLength)

            fileName = getFilename(datFiles[trajectoryNo], moleculeThiosNo)

            # plotMolecule(thiosOnly, segments, fileName)
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
        # segmentPairs = findSegmentPairs(segmentsMaster, atomMaster, maximumHoppingDistance)



        segmentPairs = [[1,2], [1,64], [1,333], [1,357], [1,356]]

        for segmentPair in segmentPairs:
            ORCAInput(segmentsMaster, segmentPair).returnAnswer()
    




        # Need to output these segments in a form that the CT program can understand
