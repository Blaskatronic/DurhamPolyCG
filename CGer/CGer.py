import os
import numpy as np
import pylab as P
import random as R
import csv
import mpl_toolkits.mplot3d.axes3d as p3
import time as T
import copy
import sys
import argparse

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
            if (datFiles[datNumber] in fileName) and ('.lammpstrj' in fileName):
                lammpsFiles.append(fileName)
                addedLammps = 1
                break
        if addedLammps == 0:
            popList.append(datNumber)
    popList = sorted(popList, reverse=True)
    for popIndex in popList:
        datFiles.pop(popIndex)
    
    return datFiles, lammpsFiles




class monomers:
    def __init__(self, originalData, bondDetails):
        self.originalData = originalData
        self.bondDetails = bondDetails
        self.originalBondDetails = copy.deepcopy(bondDetails)
        self.monomers = []
        self.organisedMonomers = []
        self.CGPointIndex = 0

    def splitMonomers(self):
        # Let's use the bond information too.
        atomDict = {'C':[], 'S':[], 'H':[]}
        totalAtoms = 0
        # Pick a sulphur to start with
        # Add both linked carbons
        # Add 2 of the three linked carbons, checking that each is not linked to a sulphur. If it is, it's the other monomer so discard
        # From those 2 carbons, add the next linked carbon (should be shared)
        # 6x (add the next linked carbon)
        # Add all hydrogens linked to all of the atoms in the monomer

        thiopheneTest = []
        alk1Test = []
        alk2Test = []
        
      
        for atom in self.originalData:
            atomDict[atom[0][0]].append(atom[1])
            totalAtoms += 1

        self.atomDictOriginal = copy.deepcopy(atomDict)

        while len(atomDict['S']) > 0:
            self.monomers.append([])
            # Select a sulphur to start with
            sulphurIndex = atomDict['S'][0]
            atomDict['S'].remove(sulphurIndex)
            self.monomers[-1].append(sulphurIndex)
            totalAtoms -= 1

            # Find the index of atoms that are linked to that sulphur
            linkedAtoms = self.findLinkedAtoms(sulphurIndex)

            for atom in linkedAtoms:
                # Iterate over all of the linked atoms sequentially
                dontAdd = 0
                # print "Checking atom", str(atom)+"..."
                # print "(linkedAtoms now reads", str(linkedAtoms)+")."
                # Find next-nearest-neighbour atoms (to check we aren't adjacent to a sulphur, which means
                # we've stepped into the next monomer)
                secondLinked = self.findLinkedAtoms(atom)
                for nextNearestAtom in secondLinked:
                    if nextNearestAtom in atomDict['S']:
                        # Don't want to add this carbon because it is from the next monomer!
                        # print "Atom is near a sulphur! Don't add! (however we just deleted the bond so we need to re-add it again)"
                        # (However, findLinkedAtoms just deleted that bond, so we need to re-add it again so we can
                        # add this carbon to its home monomer later)
                        for removedBond in self.removedBonds:
                            # print "BONDS TO ADD:", removedBond
                            self.bondDetails.append(removedBond)
                        dontAdd = 1
                if dontAdd == 0:
                    # If we're allowed to add the linked atom, add it to the right monomer, delete it from the list of
                    # untreated monomers and append the list of linkedAtoms with ITS nearest-neighbour atoms.
                    self.monomers[-1].append(atom)
                    try: 
                        atomDict['C'].remove(atom)
                    except: 
                        atomDict['H'].remove(atom)
                    totalAtoms -= 1
                    # print "Atom added successfully:"
                    # print "self.monomers =", self.monomers
                    # print "atomDict =", atomDict
                    # print "totalAtoms =", totalAtoms
                    # print "\n"
                    for nextNearestAtom in secondLinked:
                        if nextNearestAtom not in linkedAtoms:
                            linkedAtoms.append(nextNearestAtom)

        monomerDetails, monomerCGDetails, thioList, alk1List, alk2List = self.findFunctionalGroups()

        return self.monomers, self.atomDictOriginal, monomerDetails, monomerCGDetails, thioList, alk1List, alk2List



    def findLinkedAtoms(self, dummyIndex):
        attached = []
        bondsToRemove = []
        self.removedBonds = []
        for bond in self.bondDetails:
#            print bond
            if dummyIndex in bond:
                for atom in bond:
                    if atom != dummyIndex:
                        attached.append(atom)
                bondsToRemove.append(bond)
        for bond in bondsToRemove:
            # print "BONDS TO REMOVE:", bond
            self.removedBonds.append(bond)
            self.bondDetails.remove(bond)
        # print "Attached to", str(dummyIndex)+":", attached
        return attached


    def findLinkedHydrogens(self, dummyIndex, atomData):
        attached = []
        attachedDetails = []
        for bond in self.originalBondDetails:
            if dummyIndex in bond:
                for atom in bond:
                    if atom != dummyIndex:
                        if atom in self.atomDictOriginal['H']:
                            attached.append(atom)
        # Now have indices of linked hydrogens....can we get the atom data?
        for hydrogenIndex in attached:
            for atom in atomData:
                if atom[1] == hydrogenIndex:
                    attachedDetails.append(atom)
                    break
                
        return attachedDetails
        

    
    def findFunctionalGroups(self):
        thioList = []
        alk1List = []
        alk2List = []

        # Need a list of monomers including all of the atomData not just the indices
        monomerDetails = []
        monomerCGDetails = []

    #    print atomData


        for monomer in self.monomers:
    #        print "NEW MONOMER"
            monomerDetails.append([])
            monomerCGDetails.append([])
            thioList.append([])
            alk1List.append([])
            alk2List.append([])
            for index in monomer:
    #            print "Examining Index", index
                for atom in self.originalData:
    #                print "Examining atom", atom[1]
                    if atom[1] == index:
    #                    print "Indices match!"
                        monomerDetails[-1].append(atom)
                        break

    #        print monomerDetails[-1]

            thio = self.findThiophene(monomerDetails[-1])
            alk1, alk2 = self.findAlkyl(monomerDetails[-1])
            thioList[-1] = thio
            alk1List[-1] = alk1
            alk2List[-1] = alk2
            monomerCGDetails[-1] = [thio, alk1, alk2]

        return monomerDetails, monomerCGDetails, thioList, alk1List, alk2List



    def findThiophene(self, atomData):
        # Monomer comes into here to be split
        # This needed hacking - including both hydrogen atoms on the end monomers breaks LAMMPS because we have
        # 4 separate masses even though we only want 3 CG site species. As such, only add a maximum of one hydrogen.
        self.CGPointIndex += 1
        thiopheneRing = []
        hydrogensAdded = 0
        for atom in range(len(atomData)):
            if (atomData[atom][0][0] == "S") or (atomData[atom][0] == "C_2"):
                thiopheneRing.append(atomData[atom])
                hydrogens = self.findLinkedHydrogens(atomData[atom][1], atomData)
                if hydrogensAdded < 1:
                    for hydrogenAtom in hydrogens:
                        thiopheneRing.append(hydrogenAtom)
                        hydrogensAdded += 1

        thiopheneCoords, thiopheneMass = self.determineCOM(thiopheneRing)

        # Dummy index, will need changing
        thiopheneOutput = ["Thio", self.CGPointIndex, 40, thiopheneCoords, thiopheneMass]
        return thiopheneOutput


    def findAlkyl(self, atomData):
        alkyl1Carbons = []
        alkyl2Carbons = []
        C3Atoms = []
        C3AtomSep = []
        for atom in range(len(atomData)):
            if (atomData[atom][0][0] == "S"):
                sulphurPosn = atomData[atom][3]
                break

        for atom in range(len(atomData)):
            if (atomData[atom][0] == "C_3"):
                C3Atoms.append(atomData[atom])
                C3AtomSep.append(findSeparation(atomData[atom][3], sulphurPosn))


    #    print C3Atoms
    #    print C3AtomSep

        # Now have a list of the C3 atoms and a list of their separation from the Sulphur
        C3AtomSep, C3Atoms = parallelSort(C3AtomSep, C3Atoms)
        # Sorted lists. Logically, the 3 carbons closest to the sulphur are alkyl1, the 3 furthest are alkyl2

        alkyl1Details = C3Atoms[:3]
        alkyl2Details = C3Atoms[3:]

        for carbon in alkyl1Details:
            hydrogens = self.findLinkedHydrogens(carbon[1], atomData)
            for hydrogenAtom in hydrogens:
                alkyl1Details.append(hydrogenAtom)

        for carbon in alkyl2Details:
            hydrogens = self.findLinkedHydrogens(carbon[1], atomData)
            for hydrogenAtom in hydrogens:
                alkyl2Details.append(hydrogenAtom)



        alkyl1Coords, alkyl1Mass = self.determineCOM(alkyl1Details)
        alkyl2Coords, alkyl2Mass = self.determineCOM(alkyl2Details)

        # Dummy index again
        self.CGPointIndex += 1
        alkyl1 = ["Alk1", self.CGPointIndex, 24, alkyl1Coords, alkyl1Mass]
        self.CGPointIndex += 1
        alkyl2 = ["Alk2", self.CGPointIndex, 24, alkyl2Coords, alkyl2Mass]

        return alkyl1,alkyl2


    def determineCOM(self, atoms):
        massWeightedX = 0.
        massWeightedY = 0.
        massWeightedZ = 0.
        totalMass = 0.
        for atom in atoms:
            massWeightedX += atom[3][0]*atom[4]
            massWeightedY += atom[3][1]*atom[4]
            massWeightedZ += atom[3][2]*atom[4]
            totalMass += atom[4]
        return [massWeightedX/float(totalMass), massWeightedY/float(totalMass), massWeightedZ/float(totalMass)], totalMass


def loadData(filename):
    originalFile = open('./'+str(filename), 'r')
    configData = originalFile.readlines()
    originalFile.close()
    # Remove the header
    configData = configData[5:]
    # Number of atoms (assuming that only positions given!)
    numberOfAtoms = len(configData)/2
    atomDetails = []

    for i in range(numberOfAtoms):
        nameLine = treatLine(configData[2*i], 'name')
        positionLine = treatLine(configData[(2*i)+1], 'posn')
        nameLine.append(positionLine)

        # Get atom mass
        atomMass = getMass(nameLine)
        nameLine.append(atomMass)

        atomDetails.append(nameLine)

    #############################################################################
    ##    Data now in the following format (atomDetails = moleculeConfig):     ##
    ##    ["Name", index, atomicNumber, [posx, posy, posz], atomicMass]        ##
    #############################################################################

    originalField = open('./'+str(filename)+'_FIELD', 'r')
    fieldData = originalField.readlines()
    originalField.close()
    numberOfBonds = int(fieldData[6+numberOfAtoms][5:-1])
    bondDetails = []
    bondData = fieldData[7+numberOfAtoms:7+numberOfAtoms+numberOfBonds]

    for i in range(numberOfBonds):
        bond = treatLine(bondData[i], 'bond')
        bondDetails.append(bond)
    return atomDetails, bondDetails

def getMass(data):
    atomName = data[0]
    if atomName[0] == 'C':
        return 12.01115
    elif atomName[0] == 'S':
        return 32.064
    elif atomName[0] == 'H':
        return 1.00797


def treatLine(dummyLine, typeFlag):
    # Remove the spaces between the entries (Fortran formatting gubbins)
    dummyLine = dummyLine.split(' ')
    removes = []
    for position in range(len(dummyLine)):
        if len(dummyLine[position]) == 0:
            removes.append(position)
    removes = sorted(removes, reverse = True)
    for position in removes:
        dummyLine.pop(position)
    # Remove the "\n" from the end of the final entry
    dummyLine[-1] = dummyLine[-1][:-1]
    if typeFlag == 'name':
        dummyLine[1] = int(dummyLine[1])
        dummyLine[2] = int(dummyLine[2])
    elif typeFlag == 'posn':
        for i in range(len(dummyLine)):
            dummyLine[i] = float(dummyLine[i])
    elif typeFlag == 'bond':
        dummyLine = [int(dummyLine[1]), int(dummyLine[2])]
    return dummyLine


def foldedChainDict():
    # Examines a lookup file which contains the timestep number for the best approximation of a folded chain for that length
    timestepsData = {}
    timestepsHandle = open('./lookUp/foldedTimesteps.txt')
    timestepsDataRaw = timestepsHandle.readlines()
    timestepsHandle.close()
    for line in timestepsDataRaw[:-1]:
        testLine = []
        tempLine = line.split(' ')
        for element in tempLine:
            if len(element) != 0:
                testLine.append(int(element))
        timestepsData[testLine[0]] = testLine[1]
    return timestepsData



def createChain(monomersInThisChain, totalNumberOfAtoms, occupationMatrix, systemDimensionX, systemDimensionY, systemDimensionZ, chainParameters, foldedChainDictionary, rotation=False):
    # masterCGDetails needs to have the form:
    # [(outer) [(monomer) [(atom) TypeName, AtomNumber, AtomicNumber (40 for thio, 24 for alk), [posx, posy, posz], mass (81.11657 for thio, 42.08127 for Alk1, 43.08924 for Alk2) ] ] ]
    masterCGDetails = []
    monomersFileName = str(int(monomersInThisChain))
    while (len(monomersFileName) < 4):
        monomersFileName = '0'+monomersFileName

    # Load the correct lammpstrj file first:
    filename = './lookUp/'+monomersFileName+'P3HT_1.dat.lammpstrj'
    timestepData, timestepNos, simVolData = loadTrajectory(filename)

    if chainParameters[0] == True:
        # Straight chains only
        trajData = timestepData[0]
    elif chainParameters[1] == True:
        # Contorted chains only
        trajData = timestepData[1]
    elif chainParameters[2] == True:
#        timestepNumber = foldedChainDictionary[int(monomersInThisChain)]
#        trajData = timestepData[timestepNumber]
        # Folded chains only
        trajData = timestepData[4]
        pass
    elif chainParameters[3] == True:
        # Coiled chains only
        trajData = timestepData[-1]
    else:
        # Mixed chain species
        timestepNumber = R.randrange(1,11,1)
        trajData = timestepData[timestepNumber]
    # Initial Configuration = timestepData[-1]


    # HERE IS WHERE WE PICK THE TIMESTEPS.
    # Each trajectory file contains 3 distinct phases.
    # 1) The chain is 'straight' and has not folded (but is contorted as per the force field)
    # 2) The chain has 'folded' into one or more hairpins which have `stuck' together, but the chain still has linearity
    # 3) The chain has 'coiled' with multiple folds breaking linearity and resulting in an amorphous blob.

    # 'Straight' chains are always the first timestep: timestepData[0]
    # 'Contorted' chains are always the second timestep: timestepData[1]
    # 'Coiled' chains are always the final timestep: timestepData[-1]
    # 'Folded' chains are more complicated and require hard-coding in for each monomer length
    


    trajData.sort(key=lambda x: x[0]) 
    numberOfAtomsInThisChain = len(trajData)

    for atomID in np.arange(1, numberOfAtomsInThisChain, 3):
        # Of the form: [AtomID, AtomType, MoleculeID, PosX, PosY, PosZ, 0, 0, 0]
        monomerDetails = []
        totalNumberOfAtoms += 1
        monomerDetails.append(['Thio', totalNumberOfAtoms, 40, [trajData[atomID-1][3], trajData[atomID-1][4], trajData[atomID-1][5]], 81.11657])
        totalNumberOfAtoms += 1
        monomerDetails.append(['Alk1', totalNumberOfAtoms, 24, [trajData[atomID][3], trajData[atomID][4], trajData[atomID][5]], 42.08127])
        totalNumberOfAtoms += 1
        monomerDetails.append(['Alk2', totalNumberOfAtoms, 24, [trajData[atomID+1][3], trajData[atomID+1][4], trajData[atomID+1][5]], 43.08924])
        masterCGDetails.append(monomerDetails)

    # Apply rotation if required
    # And find the full extent of the randomly placed molecule to ensure that no other chains cross or intersect

    if rotation == True:
        thetaX = R.uniform(0,2*np.pi)
        thetaY = R.uniform(0,2*np.pi)
        thetaZ = R.uniform(0,2*np.pi)
    else:
        thetaX = 0
        thetaY = 0
        thetaZ = 0

    Rx = np.matrix([[1, 0, 0], [0, np.cos(thetaX), np.sin(thetaX)], [0, -np.sin(thetaX), np.cos(thetaX)]])
    Ry = np.matrix([[np.cos(thetaY), 0, -np.sin(thetaY)], [0, 1, 0], [np.sin(thetaY), 0, np.cos(thetaY)]])
    Rz = np.matrix([[np.cos(thetaZ), np.sin(thetaZ), 0], [-np.sin(thetaZ), np.cos(thetaZ), 0], [0, 0, 1]])

    rotationMatrix = Rz*Ry*Rx
    xMin = systemDimensionX
    yMin = systemDimensionY
    zMin = systemDimensionZ
    xMax = 0
    yMax = 0
    zMax = 0

    occupationCoordinates = []

    for monomerNo in range(len(masterCGDetails)):
        for CGAtomNo in range(len(masterCGDetails[monomerNo])):
            rotatedCoords = rotationMatrix*np.transpose(np.matrix(masterCGDetails[monomerNo][CGAtomNo][3]))
            newCoords = []
            roundedCoords = []
            for element in rotatedCoords:
                newCoords.append(float(element))
                roundedCoords.append(int(np.round(element)))
            if newCoords[0] < xMin:
                xMin = int(newCoords[0])
            if newCoords[0] > xMax:
                xMax = int(np.ceil(newCoords[0]))
            if newCoords[1] < yMin:
                yMin = int(newCoords[1])
            if newCoords[1] > yMax:
                yMax = int(np.ceil(newCoords[1]))
            if newCoords[2] < zMin:
                zMin = int(newCoords[2])
            if newCoords[2] > zMax:
                zMax = int(np.ceil(newCoords[2]))
            masterCGDetails[monomerNo][CGAtomNo][3] = newCoords
            occupationCoordinates.append(roundedCoords)

    # Now find a location for the COM, comparing it to the occupation matrix
    fullFlag = 0
    tooLongFlag = 0
    attempts = 0
    print "Attempting to place chain..."
    while True: # Only try to fit the polymer in 50 times, otherwise there is no space
        attempts += 1
        if np.remainder(attempts, 10) == 0:
            print "Attempt =", attempts
        if attempts == 51:
            print "Cannot fit this polymer in."
            fullFlag = 1
            break
        COMPosX = 0
        COMPosY = 0
        COMPosZ = 0
        if ((xMax-xMin) >= 2*systemDimensionX) or ((yMax-yMin) >= 2*systemDimensionY) or ((zMax-zMin) >= 2*systemDimensionZ):
            print "Simulation volume side length XYZ =", 2*systemDimensionX, 2*systemDimensionY, 2*systemDimensionZ
            print "xMax =", xMax, "xMin =", xMin, "X extent =", xMax-xMin
            print "yMax =", yMax, "yMin =", yMin, "Y extent =", yMax-yMin
            print "zMax =", zMax, "zMin =", zMin, "Z extent =", zMax-zMin
            tooLongFlag = 1
            break

        # Need to make sure the whole chain is within the simulation volume
#        print "Finding X-position within the simulation volume for the chain..."
        while True:
#            COMPosX = R.uniform(-1, 1)*(systemDimension)
            COMPosX = R.uniform(-(systemDimensionX - abs(xMin)), (systemDimensionX - abs(xMax)))
            if (abs(COMPosX + xMax) < systemDimensionX) and (abs(COMPosX + xMin) < systemDimensionX):
                break
            else:
                print "System Dimension =", systemDimensionX
                print "XMin =", xMin
                print "XMax =", xMax
                print "COMPosX =", COMPosX
                raise SystemError('INCORRECT')
        while True:
#            COMPosY = R.uniform(-1, 1)*(systemDimension)
            COMPosY = R.uniform(-(systemDimensionY - abs(yMin)), (systemDimensionY - abs(yMax)))
            if (abs(COMPosY + yMax) < systemDimensionY) and (abs(COMPosY + yMin) < systemDimensionY):
                break
            else:
                print "System Dimension =", systemDimensionY
                print "YMin =", yMin
                print "YMax =", yMax
                print "COMPosY =", COMPosY
                raise SystemError('INCORRECT')
        while True:
#            COMPosZ = R.uniform(-1, 1)*(systemDimension)
            COMPosZ = R.uniform(-(systemDimensionZ - abs(zMin)), (systemDimensionZ - abs(zMax)))
            if (abs(COMPosZ + zMax) < systemDimensionZ) and (abs(COMPosZ + zMin) < systemDimensionZ):
                break
            else:
                print "System Dimension =", systemDimensionZ
                print "ZMin =", zMin
                print "ZMax =", zMax
                print "COMPosZ =", COMPosZ
                raise SystemError('INCORRECT')



        posn = [COMPosX, COMPosY, COMPosZ]

        positionOccupied = 0

        for monomerNo in range(len(masterCGDetails)):
            for CGAtomNo in range(len(masterCGDetails[monomerNo])):
                xCheck = int(masterCGDetails[monomerNo][CGAtomNo][3][0] + COMPosX) + int(systemDimensionX)
                yCheck = int(systemDimensionY) - int(masterCGDetails[monomerNo][CGAtomNo][3][1] + COMPosY)
                zCheck = int(systemDimensionZ) - int(masterCGDetails[monomerNo][CGAtomNo][3][2] + COMPosZ)

                if xCheck >= 2*int(systemDimensionX):
                    xCheck = (2*systemDimensionX)-1
                if yCheck >= 2*int(systemDimensionY):
                    yCheck = (2*systemDimensionY)-1
                if zCheck >= 2*int(systemDimensionZ):
                    zCheck = (2*systemDimensionZ)-1
                if occupationMatrix.readValue([xCheck, yCheck, zCheck]) == 1:
                    positionOccupied = 1
                    break
            if positionOccupied == 1:
                break

        if positionOccupied == 1:
            print "Position occupied."
            continue

        else:
#            print "Increase Occupation Matrix!"
#            print "Sum before =", np.sum(occupationMatrix)


            # print occupationMatrix[int(COMPosX)+xMin+int(systemDimension):int(np.ceil(COMPosX))+xMax+int(systemDimension), int(COMPosY)+yMin+int(systemDimension):int(np.ceil(COMPosY))+yMax+int(systemDimension), int(COMPosZ)+zMin+int(systemDimension):int(np.ceil(COMPosZ))+zMax+int(systemDimension)]
            print "Moving chain to COM coordinates and updating occupation matrix..."
            for monomerNo in range(len(masterCGDetails)):
                for CGAtomNo in range(len(masterCGDetails[monomerNo])):
                    masterCGDetails[monomerNo][CGAtomNo][3][0] += COMPosX
                    masterCGDetails[monomerNo][CGAtomNo][3][1] += COMPosY
                    masterCGDetails[monomerNo][CGAtomNo][3][2] += COMPosZ

            # When updating the occupation matrix, have a 5x5x5 cube around each coordinate that is classed as "occupied"
            # REMEMBER THAT THE OCCUPATION MATRIX HAS ABSOLUTE REFERENCES (i.e. the origin is in the top left hand corner and not in the centre)
            # As such, we need to rescale all of the COM positions in order to properly account for this

            # If we assume that the origin of the matrix is at the front-top left of the occupationMatrix cube:
            # x => dimension + x
            # y => y - dimension
            # z => z - dimension

                    occupationX = int(masterCGDetails[monomerNo][CGAtomNo][3][0]) + int(systemDimensionX)
                    occupationY = int(systemDimensionY) - int(masterCGDetails[monomerNo][CGAtomNo][3][1])
                    occupationZ = int(systemDimensionZ) - int(masterCGDetails[monomerNo][CGAtomNo][3][2])
                    if occupationX >= 2*int(systemDimensionX)-10:
                        occupationX = (2*systemDimensionX)-11
                    if occupationY >= 2*int(systemDimensionY)-10:
                        occupationY = (2*systemDimensionY)-11
                    if occupationZ >= 2*int(systemDimensionZ)-10:
                        occupationZ = (2*systemDimensionZ)-11

                    if occupationX < 9:
                        occupationX = 9
                    if occupationY < 9:
                        occupationY = 9
                    if occupationZ < 9:
                        occupationZ = 9

                    for x in range(occupationX-9, occupationX+10):
                        for y in range(occupationY-9, occupationY+10):
                            for z in range(occupationZ-9, occupationZ+10):
                                occupationMatrix.insertValue([x, y, z], True)



#            print "Simulation volume is %.2f percent full..." % (np.sum(occupationMatrix)*100/(float(2*systemDimension)**3)) # DEBUGGING ONLY - THIS IS TIME CONSUMING TO CALCULATE

            # print occupationMatrix[int(COMPosX)+xMin+int(systemDimension):int(np.ceil(COMPosX))+xMax+int(systemDimension), int(COMPosY)+yMin+int(systemDimension):int(np.ceil(COMPosY))+yMax+int(systemDimension), int(COMPosZ)+zMin+int(systemDimension):int(np.ceil(COMPosZ))+zMax+int(systemDimension)]
            print "Chain placed."


            break

    
    return masterCGDetails, fullFlag, tooLongFlag, totalNumberOfAtoms, occupationMatrix


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
            if (datFiles[datNumber] in fileName) and ('.lammpstrj' in fileName):
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
    originalFile = open('./'+str(filename), 'r')
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
                if (i == 3) or (i == 4) or (i == 5):
                    tempData[i] = float(tempData[i])
                else:
                    try:
                        tempData[i] = int(tempData[i])
                    except:
                        print "Line Number", lineNumber
                        print tempData
                        raise SystemError('FAIL')
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
                    # Have to put in a check because sometimes LAMMPS outputs an exponential for very small coordinates without the actual exponent
                    # i.e. "5.0111e" which crashes the program. Something to do with the number of output characters I think, this should catch it.
                    try:
                        newData[atomNo][valueNo] = float(newData[atomNo][valueNo])
                    except ValueError:
                        newData[atomNo][valueNo] = float(0)
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


def calculateDensity(masses, traj, volume):
    # LAMMPSTRJ is written in the following format:
    # [AtomID, AtomType, MolID, x, y, z]

    # Masses is written in the following format:
    # [ [AtomType, Mass] , [AtomType2, Mass2] ... ]

    # SimVolData is written in the following format:
    # [ [ xlo, xhi ], [ ylo, yhi ], [ zlo, zhi ] ]

    totalMass = 0

    for atom in traj:
        atomType = atom[1]
        massOfAtom = masses[atomType-1][1] # g/mol
        totalMass += massOfAtom*1.660538921E-24 # g

    # Total Mass is now currently in g for the whole system volume for this timestep
    # Now calculate the volume:

    xdim = (abs(volume[0][0]) + abs(volume[0][1]))*1E-8 # cm
    ydim = (abs(volume[1][0]) + abs(volume[1][1]))*1E-8 # cm
    zdim = (abs(volume[2][0]) + abs(volume[2][1]))*1E-8 # cm

    totalVol = float(xdim*ydim*zdim) # cm^{3}

    density = totalMass/totalVol

    return density


def chainMasterFromLammpstrj(timestepData, simDims):
    # TimestepData of the form: [AtomID, AtomType, MoleculeID, PosX, PosY, PosZ, BoxX, BoxY, BoxZ]
    timestepData.sort(key = lambda x: x[0])
    numberOfAtomsInThisSnapshot = len(timestepData)
    numberOfChainsInSnapshot = 0
    totalNumberOfAtoms = 0 # running total
    chainMaster = []
    imageMaster = []
    for atom in timestepData:
        if atom[2] > numberOfChainsInSnapshot: # Molecule-No
            numberOfChainsInSnapshot = atom[2]
    for i in range(numberOfChainsInSnapshot):
        chainMaster.append([])
        imageMaster.append([])

    # minimumX = 0
    # maximumX = 0
    # minimumY = 0
    # maximumY = 0
    # minimumZ = 0
    # maximumZ = 0

    for atomID in np.arange(1, numberOfAtomsInThisSnapshot, 3):
        monomerDetails = []
        monomerImages = []
        totalNumberOfAtoms += 1
        thioX = (timestepData[atomID-1][3])#+(timestepData[atomID-1][6]*(abs(simDims[0][0])+abs(simDims[0][1])))
        thioY = (timestepData[atomID-1][4])#+(timestepData[atomID-1][7]*(abs(simDims[0][0])+abs(simDims[0][1])))
        thioZ = (timestepData[atomID-1][5])#+(timestepData[atomID-1][8]*(abs(simDims[0][0])+abs(simDims[0][1])))
        monomerDetails.append(['Thio', totalNumberOfAtoms, 40, [thioX, thioY, thioZ], 81.11657])
        monomerImages.append([timestepData[atomID-1][6], timestepData[atomID-1][7], timestepData[atomID-1][8]])
        totalNumberOfAtoms += 1
        alk1X = (timestepData[atomID][3])#+(timestepData[atomID][6]*(abs(simDims[0][0])+abs(simDims[0][1])))
        alk1Y = (timestepData[atomID][4])#+(timestepData[atomID][7]*(abs(simDims[0][0])+abs(simDims[0][1])))
        alk1Z = (timestepData[atomID][5])#+(timestepData[atomID][8]*(abs(simDims[0][0])+abs(simDims[0][1])))
        monomerDetails.append(['Alk1', totalNumberOfAtoms, 24, [alk1X, alk1Y, alk1Z], 42.08127])
        monomerImages.append([timestepData[atomID][6], timestepData[atomID][7], timestepData[atomID][8]])
        totalNumberOfAtoms += 1
        alk2X = (timestepData[atomID+1][3])#+(timestepData[atomID+1][6]*(abs(simDims[0][0])+abs(simDims[0][1])))
        alk2Y = (timestepData[atomID+1][4])#+(timestepData[atomID+1][7]*(abs(simDims[0][0])+abs(simDims[0][1])))
        alk2Z = (timestepData[atomID+1][5])#+(timestepData[atomID+1][8]*(abs(simDims[0][0])+abs(simDims[0][1])))
        monomerDetails.append(['Alk2', totalNumberOfAtoms, 24, [alk2X, alk2Y, alk2Z], 43.08924])
        monomerImages.append([timestepData[atomID+1][6], timestepData[atomID+1][7], timestepData[atomID+1][8]])

        # All atoms in this monomer will be from the same molecule, so just add them all now from the thiophene index
        moleculeIndex = timestepData[atomID-1][2]-1
        chainMaster[moleculeIndex].append(monomerDetails)
        imageMaster[moleculeIndex].append(monomerImages)

    #     if thioX > maximumX:
    #         maximumX = thioX
    #     if alk1X > maximumX:
    #         maximumX = alk1X
    #     if alk2X > maximumX:
    #         maximumX = alk2X

    #     if thioX < minimumX:
    #         minimumX = thioX
    #     if alk1X < minimumX:
    #         minimumX = alk1X
    #     if alk2X < minimumX:
    #         minimumX = alk2X

    #     if thioY > maximumY:
    #         maximumY = thioY
    #     if alk1Y > maximumY:
    #         maximumY = alk1Y
    #     if alk2Y > maximumY:
    #         maximumY = alk2Y

    #     if thioY < minimumY:
    #         minimumY = thioY
    #     if alk1Y < minimumY:
    #         minimumY = alk1Y
    #     if alk2Y < minimumY:
    #         minimumY = alk2Y

    #     if thioZ > maximumZ:
    #         maximumZ = thioZ
    #     if alk1Z > maximumZ:
    #         maximumZ = alk1Z
    #     if alk2Z > maximumZ:
    #         maximumZ = alk2Z

    #     if thioZ < minimumZ:
    #         minimumZ = thioZ
    #     if alk1Z < minimumZ:
    #         minimumZ = alk1Z
    #     if alk2Z < minimumZ:
    #         minimumZ = alk2Z

    # newSimDims = [[minimumX, maximumX], [minimumY, maximumY], [minimumZ, maximumZ]]
    # newSimDimsAbsolute = [abs(minimumX), abs(maximumX), abs(minimumY), abs(maximumY), abs(minimumZ), abs(maximumZ)]
    # newSystemDimension = max(newSimDimsAbsolute)

    return numberOfChainsInSnapshot, chainMaster, imageMaster


def duplicateData(baseP3HTCGDetails, totalChains, totalMonomersPerChain, COMPosn, totalNumberOfAtoms, rotation=False):
    masterCGDetails = []

    # How many times do we need to duplicate 2P3HT?
    duplicateNo = int(np.ceil(totalMonomersPerChain/float(len(baseP3HTCGDetails))))

    # Because we always have the same input 2P3HT, we want to translate a duplicate 2P3HT by a sensible amount
    # and then bind CG atom #1 to #4.

    # Sensible amount =  2*distance from #1-#4 (length of 2 P3HT backbone)
    # direction = along the polymer backbone

    posnOfThio1 = baseP3HTCGDetails[0][0][3]
    posnOfThio4 = baseP3HTCGDetails[1][0][3]
#    posnOfThio7 = baseP3HTCGDetails[2][0][3]
#    posnOfThio10 = baseP3HTCGDetails[3][0][3]
#    thioSeparations = []
#    thioSeparations.append(findSeparation(posnOfThio1, posnOfThio4))
#    thioSeparations.append(findSeparation(posnOfThio4, posnOfThio7))
#    thioSeparations.append(findSeparation(posnOfThio7, posnOfThio10))
#    averageThioSep = np.average(thioSeparations)
    averageThioSep = findSeparation(posnOfThio1, posnOfThio4)

    for i in range(duplicateNo):
        newMonomerSet = copy.deepcopy(baseP3HTCGDetails)
        if i != 0:
            posnOfFinalThio = masterCGDetails[-1][0][3]
            currentP3HTChainLength = findSeparation(posnOfThio1, posnOfFinalThio)
            translationDistance = currentP3HTChainLength + averageThioSep
            translationVector = findTranslationVector(posnOfThio1, posnOfFinalThio, currentP3HTChainLength, translationDistance)

        for monomerNo in range(len(newMonomerSet)):
            for CGAtomNo in range(len(newMonomerSet[monomerNo])):
                totalNumberOfAtoms += 1
                newMonomerSet[monomerNo][CGAtomNo][1] += totalNumberOfAtoms
                if i != 0:
                    newMonomerSet[monomerNo][CGAtomNo][3] = [sum(x) for x in zip(newMonomerSet[monomerNo][CGAtomNo][3], translationVector)]
            masterCGDetails.append(newMonomerSet[monomerNo])

    # Reset the centre of mass to the randomly assigned COMPosn

    # First calculate the COM translation vector
    oldCOMPosn = COMOfChain(masterCGDetails)
    COMTranslationVectorX = COMPosn[0] - oldCOMPosn[0]
    COMTranslationVectorY = COMPosn[1] - oldCOMPosn[1]
    COMTranslationVectorZ = COMPosn[2] - oldCOMPosn[2]

    # Secondly apply a random set of rotations to change the orientation of the chain

    if rotation == True:
        thetaX = R.uniform(0,2*np.pi)
        thetaY = R.uniform(0,2*np.pi)
        thetaZ = R.uniform(0,2*np.pi)
    else:
        thetaX = 0
        thetaY = 0
        thetaZ = 0

    # print "a =", np.cos(thetaX)
    # print "b =", np.sin(thetaX)
    # print "\n"
    # print "c =", np.cos(thetaY)
    # print "d =", np.sin(thetaY)
    # print "\n"
    # print "e =", np.cos(thetaZ)
    # print "f =", np.sin(thetaZ)
    # print "\n"


    Rx = np.matrix([[1, 0, 0], [0, np.cos(thetaX), np.sin(thetaX)], [0, -np.sin(thetaX), np.cos(thetaX)]])
    Ry = np.matrix([[np.cos(thetaY), 0, -np.sin(thetaY)], [0, 1, 0], [np.sin(thetaY), 0, np.cos(thetaY)]])
    Rz = np.matrix([[np.cos(thetaZ), np.sin(thetaZ), 0], [-np.sin(thetaZ), np.cos(thetaZ), 0], [0, 0, 1]])

    rotationMatrix = Rz*Ry*Rx

    for monomerNo in range(len(masterCGDetails)):
        for CGAtomNo in range(len(masterCGDetails[monomerNo])):
            masterCGDetails[monomerNo][CGAtomNo][3][0] += COMTranslationVectorX
            masterCGDetails[monomerNo][CGAtomNo][3][1] += COMTranslationVectorY
            masterCGDetails[monomerNo][CGAtomNo][3][2] += COMTranslationVectorZ
            rotatedCoords = rotationMatrix*np.transpose(np.matrix(masterCGDetails[monomerNo][CGAtomNo][3]))
            newCoords = []
            for element in rotatedCoords:
                newCoords.append(float(element))
            masterCGDetails[monomerNo][CGAtomNo][3] = newCoords
    return masterCGDetails, totalNumberOfAtoms



def COMOfChain(chain):
    massWeightedX = 0.
    massWeightedY = 0.
    massWeightedZ = 0.
    totalMass = 0.
    for monomer in chain:
        for CGAtom in monomer:
            massWeightedX += CGAtom[3][0]*CGAtom[4]
            massWeightedY += CGAtom[3][1]*CGAtom[4]
            massWeightedZ += CGAtom[3][2]*CGAtom[4]
            totalMass += CGAtom[4]
    return [massWeightedX/float(totalMass), massWeightedY/float(totalMass), massWeightedZ/float(totalMass)]


def findTranslationVector(atom1, atom2, magnitude, newMagnitude):
    xdif = (atom2[0] - atom1[0])*float(newMagnitude)/float(magnitude)
    ydif = (atom2[1] - atom1[1])*float(newMagnitude)/float(magnitude)
    zdif = (atom2[2] - atom1[2])*float(newMagnitude)/float(magnitude)
    return [xdif, ydif, zdif]

    






def writeCONFIG(filename, thio, alk1, alk2):
    originalFile = open('./'+str(filename), 'r')
    originalData = originalFile.readlines()
    originalFile.close()
    header = originalData[:5]


    dataToWrite = thio+alk1+alk2

    linesToWrite = []

    dataToWrite.sort(key=lambda x: int(x[1]))

    for CGPoint in dataToWrite:
        # First deal with name line....
        writeString = ''
        name = CGPoint[0]
        while len(name) < 15:
            name += ' '
        writeString += name
        index = str(CGPoint[1])
        while len(index) < 3:
            index = ' '+str(index)
        writeString += index
        writeString += '        '
        atomicNumber = str(CGPoint[2])
        while len(atomicNumber) < 2:
            atomicNumber = ' '+str(atomicNumber)
        writeString += atomicNumber+'\n'
        linesToWrite.append(writeString)

        # Now deal with position line....
        writeString = ''
        posx = str(CGPoint[3][0])
        while len(posx) < 18:
            posx += '0'
        if posx[0] == '-':
            posx = ' '+posx+'0'
        else:
            posx = '  '+posx
        posy = str(CGPoint[3][1])
        while len(posy) < 18:
            posy += '0'
        if posy[0] == '-':
            posy = ' '+posy+'0'
        else:
            posy = '  '+posy
        posz = str(CGPoint[3][2])
        while len(posz) < 18:
            posz += '0'
        if posz[0] == '-':
            posz = ' '+posz+'0'
        else:
            posz = '  '+posz
        writeString += posx+posy+posz+'\n'
        linesToWrite.append(writeString)


    coarseGrainedFile = open('./'+str(filename)+'_CG', 'w+')
    coarseGrainedFile.writelines(header)
    coarseGrainedFile.writelines(linesToWrite)
    coarseGrainedFile.close()

    print "Data written to: ./"+str(filename)+"_CG"

    return dataToWrite
        

    # Header will be same as before
    # Config file has foillowing format:
    # Name (15 slots), index (3 slots, align number to right), 8xSpaces, protonNumber (2 slots, align right)
    # 1x Spaces, <SIGN> posx (18 slots), 1x Spaces, <SIGN> posy (18 slots), 1x Spaces, <SIGN> posz (18 slots)
#    print header


def writeFIELD(filename, thio, alk1, alk2, dataToWrite):
    originalFile = open('./'+str(filename)+'_FIELD', 'r')
    originalData = originalFile.readlines()
    originalFile.close()
    header = originalData[:5]

    header[1] = header[1][:-2]
    header[1] += '    kcal\n' # Set units to kcal

    totalCGAtoms = len(thio)+len(alk1)+len(alk2)
    linesToWrite = []
    linesToWrite.append('ATOMS     '+str(totalCGAtoms)+'\n')

    for CGAtom in dataToWrite:
        writestring = ''
        name = CGAtom[0]
        while len(name) < 10:
            name += ' '
        writestring += name
        mass = str(CGAtom[4])
        while len(mass) < 10:
            mass += '0'
        writestring += mass+'  0.00000000\n' # This second number is the atomic charge = 0
        linesToWrite.append(writestring)

    bonds = findCGBonds(thio, alk1, alk2, dataToWrite)
    linesToWrite.append('BONDS     '+str(len(bonds))+'\n')
    
    for bond in bonds:
        writestring = ''
        ################ THIS IS WHERE WE START TO INCLUDE THE BOND PARAMETERS ETC!!!! (Separate function?)
        writestring += 'quar'
        firstCGAtom = str(bond[0])
        while len(firstCGAtom) < 5:
            firstCGAtom = ' '+firstCGAtom
        writestring += firstCGAtom
        secondCGAtom = str(bond[1])
        while len(secondCGAtom) < 5:
            secondCGAtom = ' '+secondCGAtom
        writestring += secondCGAtom
        firstParam, secondParam, thirdParam, fourthParam = findBondParams(bond, dataToWrite) # Comes out in preformatted strings
        writestring += firstParam+secondParam+thirdParam+'\n'
        linesToWrite.append(writestring)

    angles = findCGAngles(bonds)
    linesToWrite.append('ANGLES     '+str(len(angles))+'\n')

    for angle in angles:
        writestring = ''
        ############### THIS IS WHERE WE START TO INCLUDE THE BOND ANGLE PARAMETERS!!!
        writestring += 'quar' # More parameters than Quartic?? BUDDHO HELP
        firstCGAtom = str(angle[0])
        while len(firstCGAtom) < 5:
            firstCGAtom = ' '+firstCGAtom
        writestring += firstCGAtom
        secondCGAtom = str(angle[1])
        while len(secondCGAtom) < 5:
            secondCGAtom = ' '+secondCGAtom
        writestring += secondCGAtom
        thirdCGAtom = str(angle[2])
        while len(thirdCGAtom) < 5:
            thirdCGAtom = ' '+thirdCGAtom
        writestring += thirdCGAtom
        firstParam, secondParam, thirdParam, fourthParam = findAngleParams(angle, dataToWrite) # Comes out in preformatted strings
        writestring += firstParam+secondParam+thirdParam+'\n'
        linesToWrite.append(writestring)

    dihedrals = findCGAngles(angles)
    linesToWrite.append('DIHEDRALS     '+str(len(dihedrals))+'\n')

    for dihedral in dihedrals:
        writestring = ''
        ############### THIS IS WHERE WE START TO INCLUDE THE DIHEDRAL ANGLE PARAMETERS!!!
        writestring += 'opls'
        firstCGAtom = str(dihedral[0])
        while len(firstCGAtom) < 5:
            firstCGAtom = ' '+firstCGAtom
        writestring += firstCGAtom
        secondCGAtom = str(dihedral[1])
        while len(secondCGAtom) < 5:
            secondCGAtom = ' '+secondCGAtom
        writestring += secondCGAtom
        thirdCGAtom = str(dihedral[2])
        while len(thirdCGAtom) < 5:
            thirdCGAtom = ' '+thirdCGAtom
        writestring += thirdCGAtom
        fourthCGAtom = str(dihedral[3])
        while len(fourthCGAtom) < 5:
            fourthCGAtom = ' '+fourthCGAtom
        writestring += fourthCGAtom
        firstParam, secondParam, thirdParam, sixthParam, seventhParam = findDihedralParams(dihedral, dataToWrite)
        fourthParam = '  0.00000000'   # This is the 1-4 electrostatic interaction scale factor == 0
        fifthParam = '  0.00000000'    # This is the 1-4 van der Waals interaction scale factor == 0
        writestring += firstParam+secondParam+thirdParam+fourthParam+fifthParam+sixthParam+seventhParam+'\n'
        linesToWrite.append(writestring)

    linesToWrite.append('finish')


    coarseGrainedFile = open('./'+str(filename)+'_FIELD_CG', 'w+')
    coarseGrainedFile.writelines(header)
    coarseGrainedFile.writelines(linesToWrite)
    coarseGrainedFile.close()




def writeLammpsData(filename, chainMaster, imageMaster, systemDimensionX, systemDimensionY, systemDimensionZ, Mn, Mw, PolyDisp, density):
#    print chainMaster
    bondsMaster = []
    anglesMaster = []
    dihedralsMaster = []
    impropersMaster = []
    moleculeNo = 0
    for molecule in chainMaster:
        moleculeNo += 1
        print "---=== Writing molecule", moleculeNo, "===---"
        thio = []
        alk1 = []
        alk2 = []
        for monomer in molecule:
            for CGAtom in monomer:
                if CGAtom[0] == 'Thio':
                    thio.append(CGAtom)
                elif CGAtom[0] == 'Alk1':
                    alk1.append(CGAtom)
                elif CGAtom[0] == 'Alk2':
                    alk2.append(CGAtom)
                else:
                    raise SystemError('Incorrect Atom Type')
        dataToWrite = thio+alk1+alk2
        dataToWrite.sort(key = lambda x: int(x[1])) # This is nearly the same as "molecule", but not sorted by
                                                    # monomers. Required for dihedrals (span multiple monomers)
        # Make header
#        print findCGBonds(thio, alk1, alk2, dataToWrite)
        print "Finding bonds...\r"
        sys.stdout.flush()
        bonds = findCGBonds(molecule)
        print "Finding angles...\r"
        sys.stdout.flush()
        angles = findCGAngles(bonds)
        print "Finding dihedrals...\r"
        sys.stdout.flush()
        dihedrals = findCGDihedrals(bonds, angles)
        dihedrals, removedDihedrals = removeImpropers(dihedrals, dataToWrite)
        print "Finding impropers...\r"
        sys.stdout.flush()
        impropers = findCGImpropers(thio, alk1, bonds, angles)

#        impropers = impropers+removedDihedrals

        bondsMaster += bonds
        anglesMaster += angles
        dihedralsMaster += dihedrals
        impropersMaster += impropers



    # print "Bonds =", len(bondsMaster), bondsMaster
    # print "Angles =", len(anglesMaster), anglesMaster
    # print "Dihedrals =", len(dihedralsMaster), dihedralsMaster
    # print "Impropers =", len(impropersMaster), impropersMaster

    dataToWrite = []
    for molecule in chainMaster:
        for monomer in molecule:
            for atom in monomer:
                dataToWrite.append(atom)

    print "All done!"
    print "Writing header..."

    header, masses = makeLammpsHeader(dataToWrite, bondsMaster, anglesMaster, dihedralsMaster, impropersMaster, systemDimensionX, systemDimensionY, systemDimensionZ, Mn, Mw, PolyDisp, density)

    print "Writing masses..."

    massLines = writeLammpsMasses(masses)
    print "Writing atoms..."
    atomLines = writeLammpsAtoms(chainMaster, imageMaster, masses)
    print "Writing bonds..."
    bondLines = writeLammpsBonds(bondsMaster, dataToWrite)
    print "Writing angles..."
    angleLines = writeLammpsAngles(anglesMaster, dataToWrite)
    print "Writing dihedrals..."
    dihedralLines = writeLammpsDihedrals(dihedralsMaster, dataToWrite)
    print "Writing impropers..."
    improperLines = writeLammpsImpropers(impropersMaster, dataToWrite)

    print "Saving file..."

    dataFile = open('./'+str(filename)+'.dat', 'w+')
    dataFile.writelines(header)
    dataFile.writelines(massLines)
    dataFile.writelines(atomLines)
    dataFile.writelines(bondLines)
    dataFile.writelines(angleLines)
    dataFile.writelines(dihedralLines)
    dataFile.writelines(improperLines)
    dataFile.close()





def writeLammpsImpropers(impropers, dataToWrite):
    improperLines = []
    improperLines.append('\nImpropers\n\n')
    for improperNo in range(len(impropers)):
        print "\rWriting improper", improperNo, "of", str(len(impropers))+"...",
        improperID = str(improperNo + 1)
        improperType = '1' # Only one type of defined dihedral
        atom1 = str(impropers[improperNo][0])
        atom2 = str(impropers[improperNo][1])
        atom3 = str(impropers[improperNo][2])
        atom4 = str(impropers[improperNo][3])
        while len(improperID) < 8:
            improperID = ' '+improperID
        while len(improperType) < 8:
            improperType = ' '+improperType
        while len(atom1) < 8:
            atom1 = ' '+atom1
        while len(atom2) < 8:
            atom2 = ' '+atom2
        while len(atom3) < 8:
            atom3 = ' '+atom3
        while len(atom4) < 8:
            atom4 = ' '+atom4
        improperLines.append(improperID+improperType+atom1+atom2+atom3+atom4+'\n')
    print "\n"
    return improperLines



def removeImpropers(dihedrals, dataToWrite):
    # Note that the dihedrals of the form P1-P1-P1-P2 are not included in any form in these calculations.
    # This is because the planarity of the rings is conserved by the improper dihedral P1-P2-P1-P1, which
    # are included in the findCGImpropers routine.
    dihedralsToRemove = []
    for dihedralNo in range(len(dihedrals)):
        atomTypes = []
        for atom in dataToWrite:
            atomsConsidered = 0
            if atom[1] == dihedrals[dihedralNo][0]:
                atomTypes.append(atom[0])
                atomsConsidered += 1
            elif atom[1] == dihedrals[dihedralNo][1]:
                atomTypes.append(atom[0])
                atomsConsidered += 1
            elif atom[1] == dihedrals[dihedralNo][2]:
                atomTypes.append(atom[0])
                atomsConsidered += 1
            elif atom[1] == dihedrals[dihedralNo][3]:
                atomTypes.append(atom[0])
                atomsConsidered += 1
            if atomsConsidered == 4:
                break
        if (atomTypes.count('Thio') == 3) and ('Alk1' in atomTypes):
            # This potential is dictated by the improper dihedral, so we should ignore it in the dihedral calculations!
            dihedralsToRemove.append(dihedralNo)
    dihedralsToRemove.sort(reverse = True)
    removedDihedrals = []
    for index in dihedralsToRemove:
        removedDihedrals.append(dihedrals.pop(index))
    return dihedrals, removedDihedrals



def writeLammpsDihedrals(dihedrals, dataToWrite):
    # Dihedral style:
    # Format == dihedral-ID, dihedral-Type, atom1, atom2, atom3, atom4
    dihedralLines = []
    dihedralLines.append('\nDihedrals\n\n')
    for dihedralNo in range(len(dihedrals)):
        print "\rWriting dihedral", dihedralNo, "of", str(len(dihedrals))+"...",
        atomTypes = []
        for atom in dataToWrite:
            if atom[1] == dihedrals[dihedralNo][0]:
                atomTypes.append(atom[0])
                break
        for atom in dataToWrite:
            if atom[1] == dihedrals[dihedralNo][1]:
                atomTypes.append(atom[0])
                break
        for atom in dataToWrite:
            if atom[1] == dihedrals[dihedralNo][2]:
                atomTypes.append(atom[0])
                break
        for atom in dataToWrite:
            if atom[1] == dihedrals[dihedralNo][3]:
                atomTypes.append(atom[0])
                break
        if (atomTypes[0] == 'Thio') and (atomTypes[1] == 'Thio') and (atomTypes[2] == 'Thio') and (atomTypes[3] == 'Thio'):
            dihedralType = '1'
        elif (atomTypes[0] == 'Alk1') and (atomTypes[1] == 'Thio') and (atomTypes[2] == 'Thio') and (atomTypes[3] == 'Alk1'):
            dihedralType = '2'
        elif (atomTypes[0] == 'Thio') and (atomTypes[1] == 'Thio') and (atomTypes[2] == 'Alk1') and (atomTypes[3] == 'Alk2'):
            dihedralType = '3'
        elif (atomTypes[0] == 'Alk2') and (atomTypes[1] == 'Alk1') and (atomTypes[2] == 'Thio') and (atomTypes[3] == 'Thio'):
            dihedralType = '4'
        else:
            print "Dihedral in question:", dihedrals[dihedralNo]
            print "with type:", atomTypes
            raise SystemError('Incorrect dihedral type')
        dihedralID = str(dihedralNo + 1)
        atom1 = str(dihedrals[dihedralNo][0])
        atom2 = str(dihedrals[dihedralNo][1])
        atom3 = str(dihedrals[dihedralNo][2])
        atom4 = str(dihedrals[dihedralNo][3])
        while len(dihedralID) < 8:
            dihedralID = ' '+dihedralID
        while len(dihedralType) < 8:
            dihedralType = ' '+dihedralType
        while len(atom1) < 8:
            atom1 = ' '+atom1
        while len(atom2) < 8:
            atom2 = ' '+atom2
        while len(atom3) < 8:
            atom3 = ' '+atom3
        while len(atom4) < 8:
            atom4 = ' '+atom4
        dihedralLines.append(dihedralID+dihedralType+atom1+atom2+atom3+atom4+'\n')
    print "\n"
    return dihedralLines



def writeLammpsAngles(angles, dataToWrite):
    # Angle style:
    # Format == angle-ID, angle-Type, atom1, atom2, atom3
    angleLines = []
    angleLines.append('\nAngles\n\n')
    for angleNo in range(len(angles)):
        print "\rWriting angle", angleNo, "of", str(len(angles))+"...",
        atomTypes = []
        for atom in dataToWrite:
            if atom[1] == angles[angleNo][0]:
                atomTypes.append(atom[0])
                break
        for atom in dataToWrite:
            if atom[1] == angles[angleNo][1]:
                atomTypes.append(atom[0])
                break
        for atom in dataToWrite:
            if atom[1] == angles[angleNo][2]:
                atomTypes.append(atom[0])
                break
        if (atomTypes[0] == 'Thio') and (atomTypes[1] == 'Thio') and (atomTypes[2] == 'Thio'):
            angleType = '1'
        elif (atomTypes[0] == 'Thio') and (atomTypes[1] == 'Alk1') and (atomTypes[2] == 'Alk2'):
            angleType = '2'
        elif (atomTypes[0] == 'Thio') and (atomTypes[1] == 'Thio') and (atomTypes[2] == 'Alk1'):
            angleType = '3'
        elif (atomTypes[0] == 'Alk1') and (atomTypes[1] == 'Thio') and (atomTypes[2] == 'Thio'):
            angleType = '4'
        else:
            print "Angle in question:", angles[angleNo]
            print "with type:", atomTypes
            raise SystemError('Incorrect angle type')
        angleID = str(angleNo + 1)
        atom1 = str(angles[angleNo][0])
        atom2 = str(angles[angleNo][1])
        atom3 = str(angles[angleNo][2])
        while len(angleID) < 8:
            angleID = ' '+angleID
        while len(angleType) < 8:
            angleType = ' '+angleType
        while len(atom1) < 8:
            atom1 = ' '+atom1
        while len(atom2) < 8:
            atom2 = ' '+atom2
        while len(atom3) < 8:
            atom3 = ' '+atom3
        angleLines.append(angleID+angleType+atom1+atom2+atom3+'\n')
    print "\n"
    return angleLines



def writeLammpsBonds(bonds, dataToWrite):
    # Bond style:
    # Format == bond-ID, bond-Type, atom1, atom2
    bondLines = []
    bondLines.append('\nBonds\n\n')
    for bondNo in range(len(bonds)):
        print "\rWriting bond", bondNo, "of", str(len(bonds))+"...",
        atomTypes = []
        for atom in dataToWrite:
            atomsConsidered = 0
            if atom[1] == bonds[bondNo][0]:
                atomTypes.append(atom[0])
                atomsConsidered += 1
            elif atom[1] == bonds[bondNo][1]:
                atomTypes.append(atom[0])
                atomsConsidered += 1
            if atomsConsidered == 2:
                break
        if (atomTypes[0] == 'Thio') and (atomTypes[1] == 'Thio'):
            bondType = '1'
        elif ('Thio' in atomTypes) and ('Alk1' in atomTypes):
            bondType = '2'
        elif ('Alk1' in atomTypes) and ('Alk2' in atomTypes):
            bondType = '3'
        else:
            raise SystemError('Incorrect bond type')
        bondID = str(bondNo + 1)
        atom1 = str(bonds[bondNo][0])
        atom2 = str(bonds[bondNo][1])
        while len(bondID) < 8:
            bondID = ' '+bondID
        while len(bondType) < 8:
            bondType = ' '+bondType
        while len(atom1) < 8:
            atom1 = ' '+atom1
        while len(atom2) < 8:
            atom2 = ' '+atom2
        bondLines.append(bondID+bondType+atom1+atom2+'\n')
    print "\n"
    return bondLines
   

def writeLammpsAtoms(chainMaster, imageMaster, masses):
    atomLines = []
    atomLines.append('\nAtoms\n\n')
    for moleculeNo in range(len(chainMaster)):
        moleculeID = str(int(moleculeNo+1))
        print "\rWriting atoms for molecule", moleculeID+"...",
#        moleculeID = '1' # I DON'T THINK THAT MOLECULE-ID IN LAMMPS DOES WHAT I THINK IT DOES
        for monomerNumber in range(len(chainMaster[moleculeNo])):
            for CGAtomNo in range(len(chainMaster[moleculeNo][monomerNumber])):
                atomID = str(int(chainMaster[moleculeNo][monomerNumber][CGAtomNo][1]))
                x = str(chainMaster[moleculeNo][monomerNumber][CGAtomNo][3][0])[:9]
                if x[0] != '-':
                    x = x[:8]
                y = str(chainMaster[moleculeNo][monomerNumber][CGAtomNo][3][1])[:9]
                if y[0] != '-':
                    y = y[:8]
                z = str(chainMaster[moleculeNo][monomerNumber][CGAtomNo][3][2])[:9]
                if z[0] != '-':
                    z = z[:8]
                atomType = str(masses.index(chainMaster[moleculeNo][monomerNumber][CGAtomNo][4])+1)
                while len(atomID) < 8:
                    atomID = ' '+atomID
                while len(moleculeID) < 8:
                    moleculeID = ' '+moleculeID
                while len(atomType) < 8:
                    atomType = ' '+atomType
                while len(x) < 12:
                    x = ' '+x
                while len(y) < 12:
                    y = ' '+y
                while len(z) < 12:
                    z = ' '+z
                writeString = atomID+moleculeID+atomType+x+y+z
                if len(imageMaster) > 0:
                    ix = str(imageMaster[moleculeNo][monomerNumber][CGAtomNo][0])
                    iy = str(imageMaster[moleculeNo][monomerNumber][CGAtomNo][1])
                    iz = str(imageMaster[moleculeNo][monomerNumber][CGAtomNo][2])
                    while len(ix) < 12:
                        ix = ' '+ix
                    while len(iy) < 12:
                        iy = ' '+iy
                    while len(iz) < 12:
                        iz = ' '+iz
                    writeString += ix+iy+iz
                atomLines.append(writeString+'\n')   
    print "\n"
    return atomLines


def writeLammpsMasses(masses):
    massLines = []
    massLines.append('Masses\n')
    massLines.append('\n')
    for i in range(len(masses)):
        massNo = str(i+1)
        while len(massNo) < 3:
            massNo = ' '+massNo
        massVal = str(masses[i])
        massLines.append(massNo+' '+massVal+'\n')
    return massLines
    

def makeLammpsHeader(dataToWrite, bonds, angles, dihedrals, impropers, systemDimensionX, systemDimensionY, systemDimensionZ, Mn, Mw, PolyDisp, density):

    headerLines = []
    headerLines.append('LAMMPS Datafile for P3HT with Mn = '+str(int(np.round(Mn)))+' Da, Mw = '+str(int(np.round(Mw)))+' Da, polydispersity = '+str(PolyDisp)+' and density = '+str(density)+' g/cm^3\n')
    headerLines.append('\n')

    numberOfAtoms = lammpsHeaderFormatting(dataToWrite, 'atoms')
    headerLines.append(numberOfAtoms)

    numberOfBonds = lammpsHeaderFormatting(bonds, 'bonds')
    headerLines.append(numberOfBonds)

    numberOfAngles = lammpsHeaderFormatting(angles, 'angles')
    headerLines.append(numberOfAngles)
    
    numberOfDihedrals = lammpsHeaderFormatting(dihedrals, 'dihedrals')
    headerLines.append(numberOfDihedrals)
    
    numberOfImpropers = lammpsHeaderFormatting(impropers, 'impropers')
    headerLines.append(numberOfImpropers)

    headerLines.append('\n')

    masses = []
    for atom in dataToWrite:
        if atom[4] not in masses:
            masses.append(atom[4])
    # masses.sort() # There are 4 different masses - in ascending order they are alk1, alk2, thiophene[middle monomer], thiophene[end monomer]

    numberOfMasses = lammpsHeaderFormatting(masses, 'atom types')
    headerLines.append(numberOfMasses)

    numberOfBondTypes = '3'
    while len(numberOfBondTypes) < 12:
        numberOfBondTypes = ' '+numberOfBondTypes
    headerLines.append(numberOfBondTypes+'  bond types\n')

    numberOfAngleTypes = '4'
    while len(numberOfAngleTypes) < 12:
        numberOfAngleTypes = ' '+numberOfAngleTypes
    headerLines.append(numberOfAngleTypes+'  angle types\n')

    numberOfDihedralTypes = '4'
    while len(numberOfDihedralTypes) < 12:
        numberOfDihedralTypes = ' '+numberOfDihedralTypes
    headerLines.append(numberOfDihedralTypes+'  dihedral types\n')

    numberOfImproperTypes = '1'
    while len(numberOfImproperTypes) < 12:
        numberOfImproperTypes = ' '+numberOfImproperTypes
    headerLines.append(numberOfImproperTypes+'  improper types\n')

    headerLines.append('\n')

    # simDims = lammpsDimensionFormatting(dataToWrite)

    systemDimensionX = str(float(systemDimensionX))
    if len(systemDimensionX) < 8:
        while len(systemDimensionX) < 8:
            systemDimensionX += "0"
    else:
        systemDimensionX = systemDimensionX[:8]

    systemDimensionY = str(float(systemDimensionY))
    if len(systemDimensionY) < 8:
        while len(systemDimensionY) < 8:
            systemDimensionY += "0"
    else:
        systemDimensionY = systemDimensionY[:8]

    systemDimensionZ = str(float(systemDimensionZ))
    if len(systemDimensionZ) < 8:
        while len(systemDimensionZ) < 8:
            systemDimensionZ += "0"
    else:
        systemDimensionZ = systemDimensionZ[:8]

    headerLines.append('-'+systemDimensionX+' '+systemDimensionX+' xlo xhi\n')
    headerLines.append('-'+systemDimensionY+' '+systemDimensionY+' ylo yhi\n')
    headerLines.append('-'+systemDimensionZ+' '+systemDimensionZ+' zlo zhi\n')

    headerLines.append('\n')

    return headerLines, masses


def lammpsHeaderFormatting(dummyData, name):
    dummyLength = str(int(len(dummyData)))
    while len(dummyLength) < 12:
        dummyLength = ' '+dummyLength
    return dummyLength+'  '+str(name)+'\n'



def findBondParams(bond, dataToWrite):
    bondType = []
    for atomIndex in bond:
        for CGAtom in dataToWrite:
            if CGAtom[1] == atomIndex:
                bondType.append(CGAtom[0])
                break
    if bondType == ['Thio', 'Alk1']: # Huang P1-P2
        firstParam = '  105.308400' # Huang c2*2
        secondParam = '  4.07170000' # Huang l0
        thirdParam = '  613.761000' # Huang c3*3
        fourthParam = ' -66.3668000' # Huang c4*4
    elif bondType == ['Alk1', 'Alk2']: # Huang P2-P3
        firstParam = '  85.4244000'
        secondParam = '  3.53790000'
        thirdParam = ' -107.812200'
        fourthParam = ' -3061.13360'
    elif bondType == ['Thio', 'Thio']: # Huang P1-P1
        firstParam = '  138.415800'
        secondParam = '  3.82830000'
        thirdParam = '  299.937600'
        fourthParam = '  224.536400'
    return firstParam, secondParam, thirdParam, fourthParam


def findAngleParams(angle, dataToWrite):
    angleType = []
    for atomIndex in angle:
        for CGAtom in dataToWrite:
            if CGAtom[1] == atomIndex:
                angleType.append(CGAtom[0])
                break
    if angleType == ['Thio', 'Alk1', 'Alk2']: # Huang P1-P2-P3
        # No constant term in the DL_POLY potential to account for C0!!!! BUDDHO HELP
        firstParam = ' -31.8646000' # Huang c2*2
        secondParam = '  3.14159265' # Huang Phi0, Assuming radians!
        thirdParam = '  0.00000000' # Huang c3*3
        fourthParam = '  451.250800' # Huang c4*4
    elif angleType == ['Thio', 'Thio', 'Alk1']: # Huang P1-P1-P2
        firstParam = '  18.2030000'
        secondParam = '  2.13490070'
        thirdParam = ' -4.83810000'
        fourthParam = ' -37.6164000'
    elif angleType == ['Alk1', 'Thio', 'Thio']: # Huang P2-P1-P1
        firstParam = '  37.6312000'
        secondParam = '  1.44802987'
        thirdParam = '  31.9014000'
        fourthParam = ' -125.427200'
    elif angleType == ['Thio', 'Thio', 'Thio']: # Huang P1-P1-P1
        firstParam = ' -8.01400000'
        secondParam = '  3.14159265'
        thirdParam = '  0.00000000'
        fourthParam = '  134.018400'
    return firstParam, secondParam, thirdParam, fourthParam


def findDihedralParams(dihedral, dataToWrite):
    dihedralType = []
    for atomIndex in dihedral:
        for CGAtom in dataToWrite:
            if CGAtom[1] == atomIndex:
                dihedralType.append(CGAtom[0])
                break
    if dihedralType == ['Thio', 'Thio', 'Thio', 'Thio']: # Huang P1-P1-P1-P1
        firstParam = '  0.71120000' # Huang c0 - c1 + c2 - c3
        secondParam = '  1.07665000' # Huang 2c1 + (3/2)c3
        thirdParam = '  0.19480000' # Huang -c2
        sixthParam = '  0.00215000' # Huang (1/2)c3
        seventhParam = '  0.00000000'
    elif dihedralType == ['Alk1', 'Thio', 'Thio', 'Alk1']: # Huang P2-P1-P1-P2
        firstParam = ' -0.48840000' # Huang c0 - c1 + c2 - c3
        secondParam = '  2.65120000' # Huang 2c1 + (3/2)c3
        thirdParam = ' -0.70510000' # Huang -c2
        sixthParam = '  0.82480000' # Huang (1/2)c3
        seventhParam = '  0.00000000'
    elif dihedralType == ['Thio', 'Thio', 'Alk1', 'Alk2'] or dihedralType == ['Thio', 'Thio', 'Thio', 'Alk1']: # Huang P1-P1-P2-P3
        # Note, this is in here because there is no data for P1-P1-P1-P2 and otherwise it crashes =( BUDDHO HELP
        firstParam = '  0.55890000' # Huang c0 - c1 + c2 - c3
        secondParam = ' -0.71730000' # Huang 2c1 + (3/2)c3
        thirdParam = '  0.04750000' # Huang -c2
        sixthParam = '  0.07430000' # Huang (1/2)c3
        seventhParam = '  0.00000000'
    elif dihedralType == ['Alk2', 'Alk1', 'Thio', 'Thio']: # Huang P3-P2-P1-P1
        firstParam = '  0.29770000' # Huang c0 - c1 + c2 - c3
        secondParam = '  0.43590000' # Huang 2c1 + (3/2)c3
        thirdParam = ' -0.49530000' # Huang -c2
        sixthParam = '  0.00470000' # Huang (1/2)c3
        seventhParam = '  0.00000000'

    return firstParam, secondParam, thirdParam, sixthParam, seventhParam


def findCGDihedrals(bonds, angles):
    # Note that the dihedrals of the form P1-P1-P1-P2 are not included in any form in these calculations.
    # This is because the planarity of the rings is conserved by the improper dihedral P1-P2-P1-P1, which
    # are included in the findCGImpropers routine.
    dihedrals = []
    allAngles = angles[:]
    for angle in angles:
        allAngles.append([angle[2], angle[1], angle[0]])
    for angle in allAngles:
        linkedAtoms = []
        CGatom = angle[-1]
        for bond in bonds: # find all other atoms bonded to this one
            if CGatom in bond:
                for atoms in bond:
                    if (atoms not in angle): # Check the bonded atom isn't already in this angle
                        linkedAtoms.append(atoms)
                        break
        if len(linkedAtoms) != 0:
            for i in range(len(linkedAtoms)):
                dihedrals.append(angle+[linkedAtoms[i]])

    # Remove any duplicates
    dihedralCopy = copy.deepcopy(dihedrals)
    for dihedralNo in range(len(dihedralCopy)):
        dihedralCopy[dihedralNo].sort()
    dihedralTemp = []
    duplicatedDihedralIndices = []
    for dihedralNo in range(len(dihedralCopy)):
        if (dihedralCopy[dihedralNo] not in dihedralTemp):
            dihedralTemp.append(dihedralCopy[dihedralNo])
        else:
            duplicatedDihedralIndices.append(dihedralNo)
    duplicatedDihedralIndices.sort(reverse=True)
    for duplicateIndex in duplicatedDihedralIndices:
        dihedrals.pop(duplicateIndex)

    return dihedrals


def findCGAngles(bonds):
    # The original version of this function didn't find all of the angles, need to look at both directions of the bond
    angles = []
    allBonds = bonds[:]
    for bond in bonds:
        allBonds.append([bond[1], bond[0]])
    for bond in allBonds:
        linkedAtoms = []
        CGatom = bond[-1]
        for bond2 in bonds: # find all other atoms bonded to this one. Bond orientation doesn't matter here so don't need to use allBonds
            if CGatom in bond2:
                for atoms in bond2:
                    if (atoms not in bond): # Check the bonded atom isn't the atom (or bond) in question
                        linkedAtoms.append(atoms)
                        break
        if len(linkedAtoms) != 0:
            for i in range(len(linkedAtoms)):
                angles.append(bond+[linkedAtoms[i]])

    # Remove any duplicates
    angleCopy = copy.deepcopy(angles)
    for angleNo in range(len(angleCopy)):
        angleCopy[angleNo].sort()
    angleTemp = []
    duplicatedAngleIndices = []
    for angleNo in range(len(angleCopy)):
        if (angleCopy[angleNo] not in angleTemp):
            angleTemp.append(angleCopy[angleNo])
        else:
            duplicatedAngleIndices.append(angleNo)
    duplicatedAngleIndices.sort(reverse=True)
    for duplicateIndex in duplicatedAngleIndices:
        angles.pop(duplicateIndex)
    return angles


def findCGBonds(molecule):
    # The molecule data comes in as a list of monomers each with the 3 CG sites in it
    # The order is Thio, Alk1, Alk2
    bonds = []
    # Deal with intra-monomer bonds first:
    for monomerNo in range(len(molecule)):
        bonds.append([molecule[monomerNo][0][1], molecule[monomerNo][1][1]]) # P1-P2
        bonds.append([molecule[monomerNo][1][1], molecule[monomerNo][2][1]]) # P2-P3
        if ((monomerNo + 1) < len(molecule)):
            bonds.append([molecule[monomerNo][0][1], molecule[monomerNo+1][0][1]]) # P1-P1
    return bonds


def findCGImpropers(thioData, alk1Data, bonds, angles):
    # Improper is a line of 3 Thios [a, b, c] with [b] bonded to an alk1 [d]
    # The improper is then written in the order [b, d, a, c]
    thio = []
    alk1 = []
    for thiophene in thioData:
        thio.append(thiophene[1])
    for alkyl1 in alk1Data:
        alk1.append(alkyl1[1])
    impropers = []

    for angle in angles:
        if (angle[0] in thio) and (angle[1] in thio) and (angle[2] in thio):
            for bond in bonds:
                if angle[1] in bond:
                    if (angle[1] == bond[0]) and (bond[1] in alk1):
                        impropers.append([angle[1], bond[1], angle[0], angle[2]])
                    elif (angle[1] == bond[1]) and (bond[0] in alk1):
                        impropers.append([angle[1], bond[0], angle[0], angle[2]])
    return impropers



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

def findSFChainLength(SFDistribution):
    # The SF distribution is a cumulative probability distribution corresponding to chain lengths of [2:2:400] (200 elements)
    chainLengths = list(np.arange(2, 401, 2))
    randomSelection = R.random()
    for elementNo in range(len(SFDistribution)):
        if randomSelection < SFDistribution[elementNo]:
            monomerIndex = elementNo
            break
        else: # Just in case R.random() returns exactly 1, in this case return a chain length of 400
            monomerIndex = elementNo
    return chainLengths[monomerIndex]


def findSeparation(atom1, atom2):
    xdif = atom1[0] - atom2[0]
    ydif = atom1[1] - atom2[1]
    zdif = atom1[2] - atom2[2]
    return np.sqrt(xdif**2 + ydif**2 + zdif**2)


def parallelSort(list1, list2):
    data = zip(list1, list2)
    data.sort()
    list1, list2 = map(lambda t: list(t), zip(*data))
    return list1, list2
    
def plotGraph(monomerDetails, chainMaster, imageMaster, systemDimensionX, systemDimensionY, systemDimensionZ, volumeList, lammpsName, debug=False):
#    print atoms
    xOrig = []
    yOrig = []
    zOrig = []
    xCG = []
    yCG = []
    zCG = []
    colourList = []

    if debug == True:
        P1Coords = []
        P2Coords = []
        P3Coords = []
   
        for monomerNumber in range(len(monomerDetails)):
            xOrig.append([])
            yOrig.append([])
            zOrig.append([])
            for atom in monomerDetails[monomerNumber]:
                position = atom[3]
                xOrig[-1].append(position[0])
                yOrig[-1].append(position[1])
                zOrig[-1].append(position[2])

    if len(imageMaster) > 0:
        volExtentX = abs(volumeList[0][0]) + abs(volumeList[0][1])
        volExtentY = abs(volumeList[1][0]) + abs(volumeList[1][1])
        volExtentZ = abs(volumeList[2][0]) + abs(volumeList[2][1])

    for moleculeNumber in range(len(chainMaster)):
        for CGmonomerNumber in range(len(chainMaster[moleculeNumber])):
            xCG.append([])
            yCG.append([])
            zCG.append([])
            for CGAtomNo in range(len(chainMaster[moleculeNumber][CGmonomerNumber])):
                position = chainMaster[moleculeNumber][CGmonomerNumber][CGAtomNo][3]
                # for coordinate in position:
                #     if abs(coordinate) > volumeList[0][1]:
                #         print "THIS COORDINATE IS OUTSIDE THE SIMULATION"
                #         print "Atom =", chainMaster[moleculeNumber][CGmonomerNumber][CGAtomNo]
                #         print "Coordinates:", position
                if len(imageMaster) > 0: # Imported a wrapped lammpstrj, so add on the image numbers
                    position[0] += imageMaster[moleculeNumber][CGmonomerNumber][CGAtomNo][0]*volExtentX
                    position[1] += imageMaster[moleculeNumber][CGmonomerNumber][CGAtomNo][1]*volExtentY
                    position[2] += imageMaster[moleculeNumber][CGmonomerNumber][CGAtomNo][2]*volExtentZ
                xCG[-1].append(position[0])
                yCG[-1].append(position[1])
                zCG[-1].append(position[2])
                if chainMaster[moleculeNumber][CGmonomerNumber][CGAtomNo][0] == 'Thio':
                    colourList.append('r')
                elif chainMaster[moleculeNumber][CGmonomerNumber][CGAtomNo][0] == 'Alk1':
                    colourList.append('b')
                elif chainMaster[moleculeNumber][CGmonomerNumber][CGAtomNo][0] == 'Alk2':
                    colourList.append('g')
                else:
                    raise SystemError('Incorrect Atom Type')

                if debug == True:
                    coords = [position[0], position[1], position[2]]
                    if CGAtom[0] == 'Thio':
                        P1Coords.append(coords)
                    elif CGAtom[0] == 'Alk1':
                        P2Coords.append(coords)
                    elif CGAtom[0] == 'Alk2':
                        P3Coords.append(coords)
                    else:
                        raise SystemError('Incorrect Atom Type')

    if debug == True:
        separationX = np.arange(2.5,4.5,0.02)
        P1P1Separations = []
        P1P2Separations = []
        P2P3Separations = []
        for i in range(len(P1Coords)):
            P1P2Separations.append(findSeparation(P1Coords[i], P2Coords[i]))
            P2P3Separations.append(findSeparation(P2Coords[i], P3Coords[i]))
            if (i+1 < len(P1Coords)):
                P1P1Separations.append(findSeparation(P1Coords[i], P1Coords[i+1]))
        print "P1-P1 Separations =", P1P1Separations
        print "P1-P2 Separations =", P1P2Separations
        print "P2-P3 Separations =", P2P3Separations

        print "\nAverage P1-P1 Separation =", np.average(P1P1Separations)
        print "\nAverage P1-P2 Separation =", np.average(P1P2Separations)
        print "\nAverage P2-P3 Separation =", np.average(P2P3Separations)

        P1P1Hist = np.histogram(P1P1Separations, bins=100, range=(2.5,4.5), normed=True)
        P1P2Hist = np.histogram(P1P2Separations, bins=100, range=(2.5,4.5), normed=True)
        P2P3Hist = np.histogram(P2P3Separations, bins=100, range=(2.5,4.5), normed=True)

        print len(P1P1Hist[1])
        print len(P1P1Hist[0])

        P.figure()
        P.plot(P1P1Hist[1][:-1], P1P1Hist[0], c='r')
        P.plot(P1P2Hist[1][:-1], P1P2Hist[0], c='g')
        P.plot(P2P3Hist[1][:-1], P2P3Hist[0], c='b')
        P.show()
        
        raise SystemError("Debug Complete, stopping.")

    fig = P.figure()
    ax = p3.Axes3D(fig)
    colours = ['r','g','b','c','m','y','k']
    for i in range(len(xOrig)):
        # colourCode = i
        # while colourCode > 6:
        #     colourCode -= 7
        ax.plot3D(xOrig[i],yOrig[i],zOrig[i],c=colours[i%6])
    for i in range(len(xCG)):
        ax.scatter(xCG[i],yCG[i],zCG[i],s=20,c=colourList[i])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim3d(-systemDimensionX, systemDimensionX)
    ax.set_ylim3d(-systemDimensionY, systemDimensionY)
    ax.set_zlim3d(-systemDimensionZ, systemDimensionZ)
    # ax.set_xlim3d(-200, 200)
    # ax.set_ylim3d(-200, 200)
    # ax.set_zlim3d(-200, 200)
    P.savefig('./'+lammpsName+'.png')
    P.show()
    print "Figure saved as ./"+lammpsName+".png"


def checkCOMPosn(posToCheck, listOfPosns):
    # Checks that a position is at least 5 angrstroms away from a list of positions in at least 2 dimensions
    # (so that we don't get chains lying on top of each other which breaks the force field)
    positionNotViable = 0
    for position in listOfPosns:
        numberOfDimensionsCorrect = 0
        if (abs(posToCheck[0] - position[0]) >= 12.0):
            numberOfDimensionsCorrect += 1
        if (abs(posToCheck[1] - position[1]) >= 12.0):
            numberOfDimensionsCorrect += 1
        if (abs(posToCheck[2] - position[2]) >= 12.0):
            numberOfDimensionsCorrect += 1
        if numberOfDimensionsCorrect < 2:
            return 1
    return 0



class sparseMatrix:
    def __init__(self):
        self.elements = {}
    
    def insertValue(self, coordinates, value):
        self.elements[tuple(coordinates)] = value

    def readValue(self, coordinates):
        try:
            value = self.elements[tuple(coordinates)]
        except KeyError:
            value = False
        return value




def findIndex(string, logical):
    '''This function returns the locations of an inputted character (logical) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == logical:
            locations.append(index)
        index += 1
    return locations
    



if __name__ == '__main__':
    # Obtain Arguments from runCGer.py
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--lengths", default=None, help="The lengths of each polymer chain to be coarse-grained. Must be supplied if creating a new sample (use runCGer.py to submit).")
    parser.add_argument("-s", "--chainSpecies", default="con", help="The species of the chains to be created when making a new morphology (options are 'con', 'fol', 'coi', 'str' and 'mix'")
    parser.add_argument("-r", "--rotation", default=1, help="Select whether to rotate chains on placement, or leave them all aligned in a pseudo-crystal")
    parser.add_argument("-x", "--mode", default=0, help="Select the operation mode. These align with the options in runCGer.py: 0 = Generate sample, 1 = Load fastshrunk, 2 = Load heated, 3 = Load cooled, 4 = Create single chain, 5 = Create Mixed Sample")
    parser.add_argument("-d", "--dat", default=1, help="Decide whether to create a LAMMPS input file (.dat) with the coarse-grained coordinates. 1 is yes, 0 is no.")
    parser.add_argument("-g", "--graph", default=1, help="Decide whether to plot a 3D graph of the location of the coarse-grained sites in python. 1 is yes, 0 is no.")
    parser.add_argument("-ml", "--morphologyLoc", default=None, help="The directory that contains the morphology files to be treated (if mode = 1, 2 or 3)")
    parser.add_argument("-md", "--morphologyDat", default=None, help="The morphology dat file to be treated (if mode = 1, 2 or 3)")
    parser.add_argument("-mt", "--morphologyTraj", default=None, help="The morphology lammpstrj file to be treated (if mode = 1, 2 or 3)")
    args = parser.parse_args()


    # Determine chain species (see documentation for more details)
    straightChainsOnly = False
    contortedChainsOnly = False
    foldedChainsOnly = False
    coiledChainsOnly = False
    mixedChainSpecies = False

    species = args.chainSpecies

    if species == 'str':
        straightChainsOnly = True
    elif species == 'fol':
        foldedChainsOnly = True
    elif species == 'coi':
        coiledChainsOnly = True
    elif species == 'mix':
        mixedChainSpecies = True
    else:
        contortedChainsOnly = True
    chainParameters = [straightChainsOnly, contortedChainsOnly, foldedChainsOnly, coiledChainsOnly, mixedChainSpecies]

    # Parse Chain lengths array...
    chainLengths = []    
    if args.lengths != None:
        for chainLength in args.lengths.split():
            valuesOnly = ""
            for character in chainLength:
                try:
                    int(character)
                    valuesOnly += character
                except:
                    continue
            chainLengths.append(int(valuesOnly))
        totalNumberOfChains = len(chainLengths)
        rotation = bool(int(args.rotation))
    # ...or instead load the morphology to be resubmitted to LAMMPS
    else: 
        datNameRaw = args.morphologyDat
        trajectoryNameRaw = args.morphologyTraj
        direc = args.morphologyLoc
        datName = direc+datNameRaw
        trajectoryName = direc+trajectoryNameRaw
    
    # Set final parameters and initialise main variables
    mode = int(args.mode)
    doMakeLammps = bool(int(args.dat))
    doPlotGraph = bool(int(args.graph))


    chainMaster = []
    imageMaster = []
    volumeList = []
    totalNumberOfAtoms = 0
    foldedChainDictionary = foldedChainDict()

    ###### NOW, Depending on mode, EXECUTE CODE

    # First set up simulation volume:
    if len(chainLengths) == 1: # Single chain
        systemDimensionX = chainLengths[0]*4
        systemDimensionY = chainLengths[0]*4
        systemDimensionZ = chainLengths[0]*4
    else:
        systemDimensionX = 4000
        systemDimensionY = 4000
        systemDimensionZ = 4000

    # Set up the chain-Master files for writing (mode dependent)

    if (mode == 0) or (mode == 5): # Create a sample as normal, or a pre-mixed sample
        print "Creating occupation Matrix..."
        occupationMatrix = sparseMatrix()
        fullFlag = 0
        failedPlacements = 0
        numberOfChains = 0
        for length in chainLengths:
            while failedPlacements < 5:
                numberOfChains += 1
                print "Creating chain", str(numberOfChains), "of", str(totalNumberOfChains)+"..."
                chainMaster.append([])
                previousTotalNumberOfAtoms = totalNumberOfAtoms
                chainDetails, fullFlag, tooLongFlag, totalNumberOfAtoms, occupationMatrix = createChain(length, totalNumberOfAtoms, occupationMatrix, systemDimensionX, systemDimensionY, systemDimensionZ, chainParameters, foldedChainDictionary, rotation)
                if fullFlag == 1: # Tried 50 times to place chain of this orientation without success - therefore return to top of while loop and try a different random orientation for this chain
                    numberOfChains -= 1
                    failedPlacements += 1
                    chainMaster.pop(-1)
                    totalNumberOfAtoms = previousTotalNumberOfAtoms
                    print "Chain did not fit in the simulation volume..."
                    print "Failed placements =", failedPlacements
#                    break
                elif tooLongFlag == 1: # Chain will never fit because the simulation volume isn't big enough, move on to next chain length.
                    print "Chain number", numberOfChains, "containing", length, "monomers is too long for the simulation volume at this orientation."
                    failedPlacements += 1
                    numberOfChains -= 1
                    chainMaster.pop(-1)
                    totalNumberOfAtoms = previousTotalNumberOfAtoms
                else:
                    failedPlacements = 0
                    chainMaster[-1] = chainDetails
                    break
            if failedPlacements == 0:
                print "Chain placed successfully."
            else:
                print "This chain cannot be placed properly."
                raise SystemError('Not all chains can properly be placed. Reconsider.')



    elif (mode == 1): # Load fastshrunk to make .dat for heating
        trajectoryData, timestepNos, simVolData = loadTrajectory(trajectoryName)
        masses, atoms, bonds, angles, dihedrals, impropers = loadDat(datName)
        # Find out the timestep to run (at which the density >1.1 g cm^-3)
        # timestepToRun = -1 # Or the last one if we don't get to that density
        densityData = []
        for trajNo in range(len(trajectoryData)):
            volumeList = simVolData[trajNo]
            traj = trajectoryData[trajNo]
            density = calculateDensity(masses, traj, volumeList)
            densityData.append(density)
            print "Density =", density, "at timestep", timestepNos[trajNo], "which is trajectory number", trajNo
            if np.round(density, decimals=1) >= 1.1:
                timestepToRun = trajNo
                break
        print "Using the system when Density =", density, "at timestep", timestepNos[trajNo], "which is trajectory number", trajNo

        numberOfChains, chainMaster, imageMaster = chainMasterFromLammpstrj(trajectoryData[timestepToRun], simVolData[timestepToRun])
        systemDimensionX = volumeList[0][1]
        systemDimensionY = volumeList[1][1]
        systemDimensionZ = volumeList[2][1]


    elif (mode == 2) or (mode == 3): # Load heated to make .dat for cooling, or cooled .dat for equilibrating
        trajectoryData, timestepNos, simVolData = loadTrajectory(trajectoryName)
        masses, atoms, bonds, angles, dihedrals, impropers = loadDat(datName)
        # Always run the final timestep
        timestepToRun = -1
        volumeList = simVolData[timestepToRun]
        traj = trajectoryData[timestepToRun]
        numberOfChains, chainMaster, imageMaster = chainMasterFromLammpstrj(trajectoryData[timestepToRun], simVolData[timestepToRun])
        systemDimensionX = volumeList[0][1]
        systemDimensionY = volumeList[1][1]
        systemDimensionZ = volumeList[2][1]

    elif (mode == 4): # Create a single chain of specific length for lookup table
        print "Loading template 2P3HT to generate chain..."
        moleculeConfig, moleculeBonds = loadData("2P3HT")
        monomerData, atomDictOriginal, monomerDetails, monomerCGDetails, baseThioList, baseAlk1List, baseAlk2List = monomers(moleculeConfig, moleculeBonds).splitMonomers()
        # To get correct numbering, need to reset the base CG data (monomerDetails) to have all atom indices = 0. We're going to then recalculate these on an atom-by-atom basis
        for monomerNo in range(len(monomerCGDetails)):
            for CGAtomNo in range(len(monomerCGDetails[monomerNo])):
                monomerCGDetails[monomerNo][CGAtomNo][1] = 0
        print "Creating chain of length", str(chainLengths[0])+"..."
        numberOfChains = 1
        chainMaster.append([])
        masterCGDetails, totalNumberOfAtoms = duplicateData(monomerCGDetails, numberOfChains, chainLengths[0], [0, 0, 0], totalNumberOfAtoms)
        chainMaster[-1] = masterCGDetails
        systemDimensionX = chainLengths[0]*4
        systemDimensionY = chainLengths[0]*4
        systemDimensionZ = chainLengths[0]*4


    chainWeights = []
    chainWeights_2 = []
    for chain in chainMaster:
        chainWeight = 0
        chainWeight_2 = 0
        for monomer in chain:
            for atom in monomer:
                chainWeight += atom[4]
        chainWeights.append(chainWeight)
        chainWeights_2.append(chainWeight**2)

    totalMass = np.sum(chainWeights) # g/mol
    totalMassG = totalMass*1.660538921E-24 # g
    totalVolume = (2*(systemDimensionX))*(2*(systemDimensionY))*(2*(systemDimensionZ)) # ang
    totalVolumeCM3 = (2*(systemDimensionX)*1E-8)*(2*(systemDimensionY)*1E-8)*(2*(systemDimensionZ)*1E-8)

    Mn = totalMass/float(numberOfChains)
    Mw = np.sum(chainWeights_2)/totalMass
    density = totalMassG/totalVolumeCM3
    print "totalMassG =", totalMassG

    print "\nMn =", Mn, "Da"
    print "Mw =", Mw, "Da"
    PolyDisp = Mw/Mn

    print "Polydispersity =", PolyDisp
    print "Film Density =", density, "g/cm^{3}"
    print "\n"

    print "Total Number of Atoms =", totalNumberOfAtoms
    print "\n"

    weight = str(int(np.round(Mw/1000.)))
    PDI = str(np.round(PolyDisp*10)/10.)[:3]

    if (mode == 0) or (mode == 5): # Normal mode == 0
        lammpsName = 'Mw_'+weight+'_PDI_'+PDI+'_'+species+'_'+str(totalNumberOfChains)
    elif mode == 1: # load fastshrunk
        lammpsName = datNameRaw[:-4]+'_Heating'
    elif mode == 2: # load heated
        # Assuming name of the form Mw_XX_PDI_X.X_CHAINTYPE_CHAINQUANTITY_THERMALTYPE_TEMPERATURE.dat
        underscoresLoc = findIndex(datNameRaw[:-4], '_')
        lammpsName = datNameRaw[:underscoresLoc[-2]]+'_Cooling'+datNameRaw[underscoresLoc[-1]:-4]
    elif mode == 3: # load cooled
        # Assuming name of the form Mw_XX_PDI_X.X_CHAINTYPE_CHAINQUANTITY_THERMALTYPE_TEMPERATURE.dat
        underscoresLoc = findIndex(datNameRaw[:-4], '_')
        lammpsName = datNameRaw[:underscoresLoc[-2]]+'_Equil'+datNameRaw[underscoresLoc[-1]:-4]
    elif mode == 4: # Single chain, mode == 4
        chainLength = str(chainLengths[0])
        while len(chainLength) < 4:
            chainLength = '0'+chainLength
        lammpsName = chainLength+'P3HT_1'

    if doMakeLammps == True:
        print "\nChains created."
        print "\nWriting LAMMPS data..."
        writeLammpsData(lammpsName, chainMaster, imageMaster, systemDimensionX, systemDimensionY, systemDimensionZ, Mn, Mw, PolyDisp, density)
        print "LAMMPS data written to", lammpsName+".dat"
        if len(chainLengths) == 1:
            print "For the molecule of length =", chainLengths[0], "monomers."

