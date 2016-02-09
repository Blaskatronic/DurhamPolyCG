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
import pickle


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
                    try:
                        newData[atomNo][valueNo] = float(newData[atomNo][valueNo])
                    except ValueError:
                        newData[atomNo][valueNo] = 0
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



def generatePeriodicity(trajectory, volume, border):
    totalAtoms = trajectory[-1][0]
    totalChains = trajectory[-1][2]
    previousChainNumber = trajectory[-1][2]
    xMinusImage = []
    xPlusImage = []
    yMinusImage = []
    yPlusImage = []
    zMinusImage = []
    zPlusImage = []

    # This used to be set up for "totalAtoms" as the atom identifier for periodic atoms, but instead I'm going to set it to just be the original atom ID.
    # This way we can identify the atoms that are close periodically more easily.

    # print "Creating -x periodic image..."
    # Create -x image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        xCoords = atom[3]-(abs(volume[0][0])+abs(volume[0][1]))
        if xCoords >= volume[0][0]-border:
            totalAtoms += 1
            imageAtom = [atom[0], atom[1], totalChains, xCoords, atom[4], atom[5], atom[6], atom[7], atom[8]]
            xMinusImage.append(imageAtom)
    # print "Done!"
    # print "Creating +x periodic image..."

    # Create +x image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        xCoords = atom[3]+(abs(volume[0][0])+abs(volume[0][1]))
        if xCoords <= volume[0][1]+border:
            totalAtoms += 1
            imageAtom = [atom[0], atom[1], totalChains, xCoords, atom[4], atom[5], atom[6], atom[7], atom[8]]
            xPlusImage.append(imageAtom)
    # print "Done!"
    # print "Concatenating x images..."
    trajectory = xMinusImage+trajectory+xPlusImage
    # print "X-axis periodic images complete!"
    # print "Creating -y periodic image..."

    # Create -y image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        yCoords = atom[4]-(abs(volume[1][0])+abs(volume[1][1]))
        if yCoords >= volume[1][0]-border:
            totalAtoms += 1
            imageAtom = [atom[0], atom[1], totalChains, atom[3], yCoords, atom[5], atom[6], atom[7], atom[8]]
            yMinusImage.append(imageAtom)
    # print "Done!"
    # print "Creating +y periodic image..."

    # Create +y image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        yCoords = atom[4]+(abs(volume[1][0])+abs(volume[1][1]))
        if yCoords <= volume[1][1]+border:
            totalAtoms += 1
            imageAtom = [atom[0], atom[1], totalChains, atom[3], yCoords, atom[5], atom[6], atom[7], atom[8]]
            yPlusImage.append(imageAtom)
    # print "Done!"
    # print "Concatenating y images..."
    trajectory = yMinusImage+trajectory+yPlusImage
    # print "X and Y-axis periodic images complete!"
    # print "Creating -z periodic image..."

    # Create -z image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        zCoords = atom[5]-(abs(volume[2][0])+abs(volume[2][1]))
        if zCoords >= volume[2][0]-border:
            totalAtoms += 1
            imageAtom = [atom[0], atom[1], totalChains, atom[3], atom[4], zCoords, atom[6], atom[7], atom[8]]
            zMinusImage.append(imageAtom)
    # print "Done!"
    # print "Creating +z periodic image..."

    # Create +z image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        zCoords = atom[5]+(abs(volume[2][0])+abs(volume[2][1]))
        if zCoords <= volume[2][1]+border:
            totalAtoms += 1
            imageAtom = [atom[0], atom[1], totalChains, atom[3], atom[4], zCoords, atom[6], atom[7], atom[8]]
            zPlusImage.append(imageAtom)
    # print "Done!"
    # print "Concatenating z images..."

    trajectory = zMinusImage+trajectory+zPlusImage
    # print "All 8 periodic images complete!"
    # print "There are now", totalChains, "chains with", totalAtoms, "atoms in the second-order system."

    # trajectory is the entire cube, which will be required in order to make subsequent periodic cubes
    # periodicCubes is just the 8 cubes surrounding the input cube. Iterate over this when trying to find the closest atom.
    periodicCubes = xMinusImage+xPlusImage+yMinusImage+yPlusImage+zMinusImage+zPlusImage

    nextOrderVolume = []
    for dimension in volume:
        nextOrderVolume.append([])
        for value in dimension:
            if value <= 0:
                nextOrderVolume[-1].append(value-border)
            else:
                nextOrderVolume[-1].append(value+border)


    # Finally, need to get rid of the bits that extend beyond the simulation because they don't respond well to the paracrystallinity calculation
    # These bits should have been accounted for using the above periodicity generation.

    popList = []
    for atomNo in range(len(trajectory)):
        if (trajectory[atomNo][3] < nextOrderVolume[0][0]) or (trajectory[atomNo][3] > nextOrderVolume[0][1]):
            popList.append(atomNo)
        elif (trajectory[atomNo][4] < nextOrderVolume[1][0]) or (trajectory[atomNo][4] > nextOrderVolume[1][1]):
            popList.append(atomNo)
        elif (trajectory[atomNo][5] < nextOrderVolume[2][0]) or (trajectory[atomNo][5] > nextOrderVolume[2][1]):
            popList.append(atomNo)

    popList.sort(reverse=True)

    for i in popList:
        trajectory.pop(i)

    # This has now removed all perioid atoms that lie outside the simulation volume.
                  

    
    return trajectory, periodicCubes, nextOrderVolume



def calculateCrystals(atoms, bonds, angles, dihedrals, impropers, trajectory, simVolData, hkl, periodicityOrder, tolerance, crystalCutOff, border, limitThiophenesExamined, thiopheneNumber, blockAdjacent):
    # Thio atom == 1
    # Thio-Thio bond == 1
    # Thio-Thio-Thio angle == 1
    # ATOM TYPE: Molecular
    # [ atomID, moleculeID, atomType, x, y, z ]
    # Other types:
    # [ ID, Type, atom1, atom2, ...]
    if periodicityOrder > 1:
        secondOrderPeriodic = False
        secondOrderPeriodicCubes = False
        secondOrderThioAtoms = []
    if periodicityOrder > 2:
        thirdOrderPeriodic = False
        thirdOrderPeriodicCubes = False
        thirdOrderThioAtoms = []

    improperAtomsList = []
    periodicThios = []

    crystalPairs = []

    # thioAtomNos = []
    # thioAtoms = []

    # for atom in atoms:
    #      if atom[2] == 1:
    #          thioAtomNos.append(atom[0])

    # thioThioBonds = []
    # for bond in bonds:
    #     if bond[1] == 1:
    #         thioThioBonds.append([bond[2], bond[3]])
    
    # thioThioThioAngle = []
    # for angle in angles:
    #     if angle[1] == 1:
    #         thioThioThioAngle.append([angle[2], angle[3], angle[4]])

    trajectory.sort(key=lambda x: x[0])

    thioAtoms = []
    for atom in trajectory:
        if atom[1] == 1:
            thioAtoms.append(atom)

    for improper in impropers:
        improperAtoms = []
        for atomNo in improper[2:]:
            improperAtoms.append(trajectory[atomNo-1])
        improperAtomsList.append(improperAtoms)

    if limitThiophenesExamined == True:
        print "Cutting down thiophene list to", thiopheneNumber, "atoms to examine..."
        while len(improperAtomsList) > thiopheneNumber:
            thioIndexToRemove = R.randint(0, len(improperAtomsList)-1)
            improperAtomsList.pop(thioIndexToRemove)

 
#### ALONG CHAIN DIRECTION
    # How best to do this? Surely if we take the direction of the chain at each thiophene (by looking at the neighbouring bond)
    # the next thiophene in the chain will always be the closest? Maybe that's not a problem?
    # We could take the start monomer to the end monomer, but a lot of the chains are extensively twisted so I'm not so sure
    # how useful that information would be.
    # According to Poelking, d_{001} is along the chain direction and varies from 20 to 40 Angstroms which is about right.
        

#### PI STACKING DIRECTION
    dhkl = []
    dhkl2 = []
    dhkl_3 = []
    dhkl2_3 = []
    thiophenesWithoutNeighbours = 0

    print "Beginning Calculations..."

    time1 = T.time()
    totalAtoms = len(improperAtomsList)
    onePercentAtoms = int(totalAtoms/100)
    atomsTreated = 0
    plotThisOne = 0

    secondOrderPeriodic, secondOrderPeriodicCubes, secondOrderVol = generatePeriodicity(trajectory, simVolData, border)

    for atom in secondOrderPeriodic:
        if atom[1] == 1:
            secondOrderThioAtoms.append(atom)

    crystalPairsInInnerVolume = 0
    crystalPairsInOuterVolume = 0



    for improper in improperAtomsList:
        if np.remainder(atomsTreated, onePercentAtoms) == 0:
            percentageComplete = int(atomsTreated/onePercentAtoms)
            time2 = T.time()
            if percentageComplete != 0 and percentageComplete <= 100:
                totalTime = (((time2-time1)*100)/percentageComplete)
                timeRemaining = totalTime - (time2-time1)
                if timeRemaining < 60:
                    scaledTimeRem = timeRemaining
                    timeUnitsRem = 'seconds.'
                elif timeRemaining < 3600:
                    scaledTimeRem = timeRemaining/60.0
                    timeUnitsRem = 'minutes.'
                elif timeRemaining < 86400:
                    scaledTimeRem = timeRemaining/3600.0
                    timeUnitsRem = 'hours.'
                else:
                    scaledTimeRem = timeRemaining/86400.0
                    timeUnitsRem = 'days.'

                timeElapsed = time2-time1
                if timeElapsed < 60:
                    scaledTimeEl = timeElapsed
                    timeUnitsEl = 'seconds.'
                elif timeElapsed < 3600:
                    scaledTimeEl = timeElapsed/60.0
                    timeUnitsEl = 'minutes.'
                elif timeElapsed < 86400:
                    scaledTimeEl = timeElapsed/3600.0
                    timeUnitsEl = 'hours.'
                else:
                    scaledTimeEl = timeElapsed/86400.0
                    timeUnitsEl = 'days.'
                print "\rCalculations %d percent completed. Remaining time estimate: %.1f %s. Elapsed time: %.1f %s." % (percentageComplete, scaledTimeRem, timeUnitsRem, scaledTimeEl, timeUnitsEl),



        if hkl == [0, 1, 0]:
            stackingVector = piStackingDirection(improper)
        elif hkl == [0, 0, 1]:
            stackingVector = alongChainDirection(improper)
        elif hkl == [1, 0, 0]:
            raise SystemError("Alkyl stacking subroutines not yet implemented")
        else:
            print "hkl =", hkl
            raise SystemError("Higher order hkl values not yet implemeneted")
        
        if blockAdjacent == True:
            bondedAtoms = []
            for bond in bonds:
                if (improper[0][0] == bond[2]):
                    bondedAtoms.append(bond[3])
                elif (improper[0][0] == bond[3]):
                    bondedAtoms.append(bond[2])
        else:
            bondedAtoms = [0]


        closestAtoms = findClosestAtoms(improper[0], secondOrderThioAtoms, stackingVector, tolerance, bondedAtoms, crystalCutOff)
        if (closestAtoms == 0):
            thiophenesWithoutNeighbours += 1
            closestAtoms = [[0, 100, 0, [0,0,0], [0,0,0]]]
        else:
            for atom in closestAtoms:
                if improper[0][0] < atom[0]:
                    if [improper[0][0], atom[0]] not in crystalPairs:
                        crystalPairs.append([improper[0][0], atom[0]])
                        atom1Posn = [improper[0][3], improper[0][4], improper[0][5]]
                        atom2Posn = atom[4]
                        averagePosition = [((atom1Posn[0] + atom2Posn[0])/2.), ((atom1Posn[1] + atom2Posn[1])/2.), ((atom1Posn[2] + atom2Posn[2])/2.)]
                        if ((averagePosition[0] < simVolData[0][1]/(2**(1/3.))) and (averagePosition[0] > simVolData[0][0]/(2**(1/3.)))):
                            if ((averagePosition[1] < simVolData[1][1]/(2**(1/3.))) and (averagePosition[1] > simVolData[1][0]/(2**(1/3.)))):
                                if ((averagePosition[2] < simVolData[2][1]/(2**(1/3.))) and (averagePosition[2] > simVolData[2][0]/(2**(1/3.)))):
                                    crystalPairsInInnerVolume += 1
                                else:
                                    crystalPairsInOuterVolume += 1
                            else:
                                crystalPairsInOuterVolume += 1
                        else:
                            crystalPairsInOuterVolume += 1


                else:
                    if [atom[0], improper[0][0]] not in crystalPairs:
                        crystalPairs.append([atom[0], improper[0][0]])
                        atom1Posn = [improper[0][3], improper[0][4], improper[0][5]]
                        atom2Posn = atom[4]
                        averagePosition = [((atom1Posn[0] + atom2Posn[0])/2.), ((atom1Posn[1] + atom2Posn[1])/2.), ((atom1Posn[2] + atom2Posn[2])/2.)]
                        if ((averagePosition[0] < simVolData[0][1]/(2**(1/3.))) and (averagePosition[0] > simVolData[0][0]/(2**(1/3.)))):
                            if ((averagePosition[1] < simVolData[1][1]/(2**(1/3.))) and (averagePosition[1] > simVolData[1][0]/(2**(1/3.)))):
                                if ((averagePosition[2] < simVolData[2][1]/(2**(1/3.))) and (averagePosition[2] > simVolData[2][0]/(2**(1/3.)))):
                                    crystalPairsInInnerVolume += 1
                                else:
                                    crystalPairsInOuterVolume += 1
                            else:
                                crystalPairsInOuterVolume += 1
                        else:
                            crystalPairsInOuterVolume += 1
                            


        dhkl.append(closestAtoms[0][1])
        dhkl2.append(closestAtoms[0][1]**2)


        atomsTreated += 1



    avdhkl = np.sum(dhkl)/atomsTreated
    avdhkl2 = np.sum(dhkl2)/atomsTreated

    g2 = (avdhkl2 - avdhkl**2)/(avdhkl**2)


    print "\n-----=====================-----"
    # print "Paracrystallinity parameter g =", np.sqrt(g2)
    # print "Proportion of thiophenes with no neighbours (calcs truncated) =", float(thiophenesWithoutNeighbours)/float(atomsTreated)
    print "Proportion of crystal pairs in inner volume compared to outer =", float(crystalPairsInInnerVolume)/float(crystalPairsInInnerVolume+crystalPairsInOuterVolume)
    print "-----=====================-----\n"

    # P.figure()
    # x = np.arange(0, int(max(dhkl)))
    # P.hist(dhkl, x)
    # P.show()

    # P.figure()
    # x = np.arange(0, int(max(dhkl2)))
    # P.hist(dhkl2, x)
    # P.show()

    return np.sqrt(g2), dhkl, crystalPairs

#    print len(atoms)
#    print len(thioAtoms)



# def findClosestAtom(thioAtom, allThioAtoms, axis, toleranceDistance): # LINEAR
#     closestAtomsList = []
#     # print "--== BEGIN OUTPUT ==--"
#     # print "ThioAtom =", thioAtom
#     # print "Axis =", axis
#     # Draw a line between the thiophene atom and the position 1 unit along the given axis
#     x1 = np.array([thioAtom[3], thioAtom[4], thioAtom[5]])
#     x2 = np.array([thioAtom[3]+axis[0], thioAtom[4]+axis[1], thioAtom[5]+axis[2]])
#     # The line x1-x2 describes the pi stacking axis from the thioAtom
#     for atom in allThioAtoms:
#         if atom == thioAtom:
#             # The atom should be present in allThioAtoms, so ignore itself
#             continue
#         x0 = np.array([atom[3], atom[4], atom[5]])
#         vec1 = x0-x1
#         vec2 = x0-x2
#         # Find the magnitude of the cross product of the two difference vectors vec1 and vec2
#         crossVec = np.cross(vec1, vec2)
#         numerator = np.linalg.norm(crossVec)
#         denominator = np.linalg.norm(x2-x1) # Should always be 1, left in to check
#         separationFromAxis = numerator/denominator
#         # print "Atom =", atom
#         # print "x1 =", x1
#         # print "x2 =", x2
#         # print "x0 =", x0
#         # print "x0 - x1 =", vec1
#         # print "x0 - x2 =", vec2
#         # print "Cross product =", crossVec
#         # print "Numerator =", numerator
#         # print "Denominator =", denominator
#         # print "Separation =", separationFromAxis
#         if separationFromAxis <= toleranceDistance:
#             # Append: [atomIndex, distanceFromThiophere, separationFromPlane]
#             closestAtomsList.append([atom[0], np.linalg.norm(vec1), separationFromAxis])
#     closestAtomsList.sort(key = lambda x: x[1])
#     if len(closestAtomsList) != 0:
#         return closestAtomsList[0]
#     else:
#         return 0


# def findClosestAtom(atom, atomList, axis):
#     # Axis describes normal vector to the plane
#     toleranceDistance = 10 # Distance a conjugated subunit can be from the plane to still be classed as in that plane (angstroms)
#     closestAtomsList = []
#     atom1Index = atom[0]
#     atom1Position = [atom[3], atom[4], atom[5]]
#     # First obtain the equation of the plane ax + by + cz + d = 0
#     # [a,b,c] = [axis[0], axis[1], axis[2]]
#     # d = - a*x0 - b*y0 - c*z0
#     # where [x0, y0, z0] is a point in the plane ([atom[0], atom[1], atom[2]])
#     # Using this naming nomenclature:
#     a = axis[0]
#     b = axis[1]
#     c = axis[2]
#     d = - (axis[0]*atom[0]) - (axis[1]*atom[1]) - (axis[2]*atom[2])
#     printing = False

#     # if atom1Index == 244:
#     #     raw_input('Ready to start printing, hit return to continue (spam inc)')
#     #     printing = True


#     for atom2 in atomList:
#         if atom == atom2:
#             continue
#         atom2Index = atom2[0]
#         atom2Position = [atom2[3], atom2[4], atom2[5]]
#         separationFromPlane = abs((a*atom2[3]) + (b*atom2[4]) + (c*atom2[5]) + d)/np.sqrt(a**2+b**2+c**2)
#         if printing == True:
#             print "Atom2 =", atom2
#             print "Separation from plane =", separationFromPlane
#         if separationFromPlane <= toleranceDistance:
#             closestAtomsList.append([atom2Index, calculateSeparation(atom1Position, atom2Position), separationFromPlane])


#     closestAtomsList.sort(key = lambda x: x[1])
#     # ClosestAtomsList is in the format: [ Atom Index, Separation, AngleBetweenSeparationVectorAndAxis, [SepVecX, SepVecY, SepVecZ] ]
#     if len(closestAtomsList) != 0:
#         return closestAtomsList[0]
#     else:
#         return 0

def getBondDict(bonds):
    bondDict = {}
    for bond in bonds:
        if bond[1] == 1:
            if bond[2] not in bondDict:
                bondDict[bond[2]] = [bond[3]]
            else:
                bondDict[bond[2]].append(bond[3])
            if bond[3] not in bondDict:
                bondDict[bond[3]] = [bond[2]]
            else:
                bondDict[bond[3]].append(bond[2])
    return bondDict


def treatCrystalPairs(rawCrystalPairs, bondDict):
    t0 = T.time()
    crystalPairs = sorted(rawCrystalPairs)
    crystalsAssorted = []
    crystals = []
    while len(crystalPairs) > 0:
        currentCrystal = []
        popList = []
        for index in crystalPairs[0]:
            currentCrystal.append(index)
        crystalPairs.pop(0)
        while True:
            addedNewIndex = False
            for crystalPairNo in range(len(crystalPairs)):
                crystalPair = crystalPairs[crystalPairNo]
                if (crystalPair[0] in currentCrystal) and (crystalPair[1] not in currentCrystal):
                    addedNewIndex = True
                    currentCrystal.append(crystalPair[1])
                    # # As a rough idea, try adding adjacent thiophenes to the crystal too, as a rough way to populate crystals with adjacent neighbours
                    # for neighbourThio in bondDict[crystalPair[1]]:
                    #     if (neighbourThio not in currentCrystal):
                    #         currentCrystal.append(neighbourThio)
                    popList.append(crystalPairNo)
                if (crystalPair[1] in currentCrystal) and (crystalPair[0] not in currentCrystal):
                    addedNewIndex = True
                    currentCrystal.append(crystalPair[0])
                    # # As above
                    # for neighbourThio in bondDict[crystalPair[0]]:
                    #     if (neighbourThio not in currentCrystal):
                    #         currentCrystal.append(neighbourThio)
                    popList.append(crystalPairNo)
            if addedNewIndex == False:
                break
        popList = sorted(popList, reverse=True)
        for elementNo in popList:
            crystalPairs.pop(elementNo)
        crystals.append(currentCrystal)
        for element in currentCrystal:
            crystalsAssorted.append(element)


    # Link together crystals that are neighbouring
    crystalsIncNeighbours = []
    while len(crystals) > 0:
        currentCrystal = copy.deepcopy(crystals[0])
        crystals.pop(0)
        crystalIDsAdded = []
        for atomID in currentCrystal:
            for neighbourThio in bondDict[atomID]:
                for remainingCrystal in crystals:
                    if neighbourThio in remainingCrystal:
                        if crystals.index(remainingCrystal) not in crystalIDsAdded:
                            crystalIDsAdded.append(crystals.index(remainingCrystal))
                            currentCrystal = currentCrystal+remainingCrystal
                            break
        crystalIDsAdded = sorted(crystalIDsAdded, reverse=True)
        for crystalID in crystalIDsAdded:
            crystals.pop(crystalID)
        crystalsIncNeighbours.append(currentCrystal)
                                                 
    t1 = T.time()                
    print "There are", len(crystalsIncNeighbours), "crystals. Elasped Time =", t1-t0, "seconds."

    return crystalsIncNeighbours, crystalsAssorted


def getBiggestCrystal(crystals):
    longestCrystalSize = 0
    longestCrystalIndex = 0
    for crystalNo in range(len(crystals)):
        if len(crystals[crystalNo]) > longestCrystalSize:
            longestCrystalSize = len(crystals[crystalNo])
            longestCrystalIndex = crystalNo
    print "The longest crystal is number", longestCrystalIndex, "with a total size of", longestCrystalSize
    return crystals[longestCrystalIndex]
        

def getCrystalsPerChain(atoms, crystals):
    crystalsPerChain = {}
    for crystal in crystals:
        chainsInThisCrystal = []
        for atomIndex in crystal:
            for atom in atoms:
                if atom[0] == atomIndex:
                    if atom[1] not in chainsInThisCrystal:
                        chainsInThisCrystal.append(atom[1])
                        break
        for chainIndex in chainsInThisCrystal:
            if chainIndex not in crystalsPerChain:
                crystalsPerChain[chainIndex] = [crystals.index(crystal)]
            else:
                crystalsPerChain[chainIndex].append(crystals.index(crystal))
        # print "Crystal Index =", crystals.index(crystal), "Chains in crystal =", chainsInThisCrystal
    return crystalsPerChain

def getChainLengths(atoms):
    chainLengths = {}
    for atom in atoms:
        if atom[1] not in chainLengths:
            chainLengths[atom[1]] = 1
        else:
            chainLengths[atom[1]] += 1
    return chainLengths

def makeCSVs(name, crystalsPerChain, crystals, chainLengths):
    crystalsPerChainFile = './outputFiles/'+name+'_crystPerChain.csv'
    csvFile = csv.writer(open(crystalsPerChainFile, 'w+'), delimiter = ',')
    for chainNo in crystalsPerChain:
        data = [chainNo]
        for crystalID in crystalsPerChain[chainNo]:
            data.append(crystalID)
        csvFile.writerow(data)
    print 'Crystals Per Chain written to', crystalsPerChainFile

    crystalsVLength = './outputFiles/'+name+'_crystalVLength.csv'
    document = csv.writer(open(crystalsVLength, 'w+'), delimiter = ',')
    for chainID in crystalsPerChain:
        data = [chainLengths[chainID]/3, len(crystalsPerChain[chainID])] # Need to divide the chain lengths by 3 to get the number of monomers, otherwise we include the alk1 and alk2 bits too.
        document.writerow(data)
    print 'Length vs Crystals data written to', crystalsVLength
    


def makeLammpsFileChainTest(name, chainNo, crystalAtomsInChain, crystals):
    currentLammpstrjHandle = open(name+'.dat.lammpstrj', 'r')
    currentLammpstrjLines = currentLammpstrjHandle.readlines()
    currentLammpstrjHandle.close()

    chainLines = []

    for i in range(9):
        chainLines.append(currentLammpstrjLines[i])

    chainAtoms = 0

    for i in range(9, len(currentLammpstrjLines)):
        testLine = currentLammpstrjLines[i]
        testLine = testLine.split(' ')
        index = int(testLine[0])
        atomType = int(testLine[1])
        moleculeNumber = int(testLine[2])
        if moleculeNumber == chainNo:
            for crystalNo in crystalAtomsInChain:
                if index in crystals[crystalNo]:
                    testLine[1] = str(4+crystalNo)
            chainAtoms += 1
            newLine = ' '.join(testLine)
            chainLines.append(newLine)
        else:
            continue

    # Now change the atom types to 4+ for each crystal so we can highlight them different colours in VMD            
        


    chainLines[3] = str(chainAtoms)+"\n"
    chainFileName = "chain_"+name+'_'+str(chainNo)+'.dat.lammpstrj'
    chainLammpsFileHandle = open('./outputFiles/'+chainFileName, 'w+')
    # Update the correct number of atoms ((number of lines to write) - 9)
    # chainAdditionLines[chainIndex][3] = str(int(len(chainAdditionLines[chainIndex])-9))+'\n'
    chainLammpsFileHandle.writelines(chainLines)
    chainLammpsFileHandle.close()




def makeNewLammpsFiles(name, crystals, biggestCrystal, chainLengthAdditions):
    currentLammpstrjHandle = open(name+'.dat.lammpstrj', 'r')
    currentLammpstrjLines = currentLammpstrjHandle.readlines()
    currentLammpstrjHandle.close()

    allCrystalsLammpstrjLines = []
    biggestCrystalLammpstrjLines = []
    sliceLammpstrjLines = []
    chainAdditionLines = {}
    if chainLengthAdditions != None:
        for addedChain in chainLengthAdditions:
            chainAdditionLines[addedChain] = []


    # Put in the header first. NOTE THAT THE NUMBER OF ATOMS WILL BE WRONG
    for i in range(9):
        allCrystalsLammpstrjLines.append(currentLammpstrjLines[i])
        biggestCrystalLammpstrjLines.append(currentLammpstrjLines[i])
        sliceLammpstrjLines.append(currentLammpstrjLines[i])
        if chainLengthAdditions != None:
            for jKey in chainAdditionLines:
                chainAdditionLines[jKey].append(currentLammpstrjLines[i])
            

    allCrystalsNumberOfAtoms = 0
    biggestCrystalNumberOfAtoms = 0
    sliceNumberOfAtoms = 0

    for i in range(9, len(currentLammpstrjLines)):
        testLine = currentLammpstrjLines[i]
        testLine = testLine.split(' ')
        index = int(testLine[0])
        moleculeNumber = int(testLine[2])
        zCoordinate = float(testLine[5])
        for crystalIndex in range(len(crystals)):
            if index in crystals[crystalIndex]:
                allCrystalsNumberOfAtoms += 1
                testLine[1] = str(crystalIndex+1)
                newLine = ' '.join(testLine)
                allCrystalsLammpstrjLines.append(newLine)
        # if index in crystals:
        #     allCrystalsNumberOfAtoms += 1
        #     allCrystalsLammpstrjLines.append(currentLammpstrjLines[i])
        if index in biggestCrystal:
            biggestCrystalNumberOfAtoms += 1
            biggestCrystalLammpstrjLines.append(currentLammpstrjLines[i])
        if (zCoordinate > -2.) and (zCoordinate < 2.):
            sliceNumberOfAtoms += 1
            sliceLammpstrjLines.append(currentLammpstrjLines[i])
        if chainLengthAdditions != None:
            if moleculeNumber in chainLengthAdditions:
                chainAdditionLines[moleculeNumber].append(currentLammpstrjLines[i])

    allCrystalsLammpstrjLines[3] = str(allCrystalsNumberOfAtoms)+"\n"
    biggestCrystalLammpstrjLines[3] = str(biggestCrystalNumberOfAtoms)+"\n"
    sliceLammpstrjLines[3] = str(sliceNumberOfAtoms)+"\n"

    allCrystalsFileName = "crys_"+name+'.dat.lammpstrj'
    
    allCrystalsLammpsFileHandle = open('./outputFiles/'+allCrystalsFileName, 'w+')
    allCrystalsLammpsFileHandle.writelines(allCrystalsLammpstrjLines)
    allCrystalsLammpsFileHandle.close()

    biggestCrystalFileName = "big_"+name+'.dat.lammpstrj'
    
    biggestCrystalLammpsFileHandle = open('./outputFiles/'+biggestCrystalFileName, 'w+')
    biggestCrystalLammpsFileHandle.writelines(biggestCrystalLammpstrjLines)
    biggestCrystalLammpsFileHandle.close()

    sliceFileName = "slice_"+name+'.dat.lammpstrj'
    
    sliceLammpsFileHandle = open('./outputFiles/'+sliceFileName, 'w+')
    sliceLammpsFileHandle.writelines(sliceLammpstrjLines)
    sliceLammpsFileHandle.close()

    if chainLengthAdditions != None:
        for chainIndex in chainAdditionLines:
            chainFileName = "chain_"+name+'_'+str(chainIndex)+'.dat.lammpstrj'
            chainLammpsFileHandle = open('./outputFiles/'+chainFileName, 'w+')
            # Update the correct number of atoms ((number of lines to write) - 9)
            chainAdditionLines[chainIndex][3] = str(int(len(chainAdditionLines[chainIndex])-9))+'\n'
            chainLammpsFileHandle.writelines(chainAdditionLines[chainIndex])
            chainLammpsFileHandle.close()
        



def findClosestAtoms(thioAtom, allThioAtoms, axis, toleranceAngle, bondedAtoms, crystalCutOff): # ANGULAR
#    toleranceAngle = np.pi/180 # Monomers subtending an angle of toleranceAngle degrees are still classed as in the same plane.
#    toleranceDistance = 2 # Distance a conjugated subunit can be from the plane to still be classed as in that plane (angstroms)
#    print "AXIS =", axis
    closestAtomsList = []
    atom1Index = thioAtom[0]
    atom1Position = [thioAtom[3], thioAtom[4], thioAtom[5]]
    for atom2 in allThioAtoms:
        if (atom2 == thioAtom) or (atom2[0] in bondedAtoms):
            continue
        atom2Index = atom2[0]
        atom2Position = [atom2[3], atom2[4], atom2[5]]
        separationVector = np.array([atom2[3] - thioAtom[3], atom2[4] - thioAtom[4], atom2[5] - thioAtom[5]])
        separationVector = list(normaliseVec(separationVector))

        # a.b = |a| |b| cos (theta)
        dotProductAngle = np.arccos(np.dot(separationVector, axis))

        # Periodicity of sin/cos means it might measure theta in the opposite direction, catch that (we want the smallest theta)
        if dotProductAngle >= (np.pi/2.):
            dotProductAngle = np.pi - dotProductAngle
        if dotProductAngle > np.pi:
            # Should be completely impossible
            print "---======---"
            print "Atom =", thioAtom
            print "Atom2 =", atom2
            print "Separation Vector =", separationVector
            print "Axis =", axis
            print "Cross product =", np.cross(separationVector, axis)
            print "Normalised Cross Product =", np.linalg.norm(np.cross(separationVector, axis))
            print "dotProductAngle =", dotProductAngle
            raise SystemError('Angle calculation is incorrect')

        if dotProductAngle <= toleranceAngle:
            # This monomer is within the axis plane therefore add it to the closestAtomsList
            separation = calculateSeparation(atom1Position, atom2Position)
            if separation <= crystalCutOff:
                # This monomer is within a crystalline proximity to the neighbour, therefore add it to the closestAtomsList
                closestAtomsList.append([atom2Index, calculateSeparation(atom1Position, atom2Position), dotProductAngle, separationVector, atom2Position])

    closestAtomsList.sort(key = lambda x: x[1])
#    print closestAtomsList

#    print "Length =", closestAtomsList

    # ClosestAtomsList is in the format: [ Atom Index, Separation, AngleBetweenSeparationVectorAndAxis, [SepVecX, SepVecY, SepVecZ] ]
    if len(closestAtomsList) != 0:
        return closestAtomsList
    else:
        return 0

def calculateSeparation(atom1, atom2):
    xdif = atom1[0] - atom2[0]
    ydif = atom1[1] - atom2[1]
    zdif = atom1[2] - atom2[2]
    return np.sqrt(xdif**2 + ydif**2 + zdif**2)

def normaliseVec(vector):
    return vector/float(np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2))

def piStackingDirection(angle):
    # Angle comes in as a 4-vector of:
    # [AtomNo, Type, Molecule, x, y, z, ix, iy, iz] of type 1-2-1-1
    # Normal direction to a plane defined by points P1, P2 and P3 =
    # (P3 - P1) vectorProduct (P2 - P1)
    # print "Atom =", atom
    vec1 = np.array([(angle[2][3] - angle[0][3]), (angle[2][4] - angle[0][4]), (angle[2][5] - angle[0][5])])
    vec2 = np.array([(angle[1][3] - angle[0][3]), (angle[1][4] - angle[0][4]), (angle[1][5] - angle[0][5])])
    normalVec = np.cross(vec1, vec2)
    normalVec = normaliseVec(normalVec)
    return normalVec


def alongChainDirection(angle):
    # Along the chain direction is the final element of the above improper (type 1-2-1-1) subtracted from the first
    # e.g. the first improper is [4, 5, 1, 7]. Take the along-chain direction to be between 4 and the next thiophene along which is 7
    chainVector = np.array([(angle[3][3] - angle[0][3]), (angle[3][4] - angle[0][4]), (angle[3][5] - angle[0][5])])
    chainVector = normaliseVec(chainVector)
    return chainVector
    


def findNearbyBoundaries(atom, simVolData):
    # Check SimVolData is ((xMin, xMax), (yMin, yMax), (zMin, zMax))
    border = 10         # The distance, beyond which an atom is said to be close to the border and we need to include extra periods
    xPosn = atom[3]
    yPosn = atom[4]
    zPosn = atom[5]
    periodsToInclude = [0, 0, 0, 0, 0, 0] # [MinusX, PlusX, MinusY, PlusY, MinusZ, PlusZ]
    if (xPosn - simVolData[0][0]) < border:
        periodsToInclude[0] = 1
    if (simVolData[0][1] - xPosn) < border:
        periodsToInclude[1] = 1
    if (yPosn - simVolData[0][0]) < border:
        periodsToInclude[2] = 1
    if (simVolData[1][1] - yPosn) < border:
        periodsToInclude[3] = 1
    if (zPosn - simVolData[2][0]) < border:
        periodsToInclude[4] = 1
    if (simVolData[2][1] - zPosn) < border:
        periodsToInclude[5] = 1
    return periodsToInclude
    


def findIndex(string, logical):
    '''This function returns the locations of an inputted character (logical) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == logical:
            locations.append(index)
        index += 1
    return locations
    
def writeCrystalHistCSV(data, name):
    '''Writes the dhkl values in a CSV file'''
    filename = './outputFiles/'+name+'_crystalHist.csv'
    document = csv.writer(open(filename, 'w+'), delimiter = ',')
    for dhkl in data:
        document.writerow([dhkl])
    print 'Crystals per chain histogram written to', filename


def plot(timestepNos, paraData, filename):
    P.figure()
    P.plot(timestepNos, paraData)
    P.xlabel('Timestep')
    P.ylabel('Paracrystallinity, g')
    P.savefig('./outputGraphs/para_'+str(filename)+'.png')


def plotHist(data, name):
    thiosUnder7 = []
    for dhkl in data:
        if dhkl < 7:
            thiosUnder7.append(dhkl)


    P.figure()
    n, bins, patches = P.hist(data, 200, normed=1)
    P.xlabel('dhkl')
    P.ylabel('Occurence')
    P.xlim(0, 100)
    P.ylim(0.08)
    P.savefig('./outputGraphs/'+name+'.png')


    n, bins, patches = P.hist(thiosUnder7, 10, normed=1)

    minVal = 99999999
    mind = 0
    for i in range(len(n)):
        if (n[i] != 0) and (n[i] < minVal):
            minVal = n[i]
            mind = bins[i+1]

    thiosInCrystal = thiosUnder7

    popList = []
    for dhklNo in range(len(thiosInCrystal)):
        if thiosInCrystal[dhklNo] > mind:
            popList.append(dhklNo)
    popList = sorted(popList, reverse=True)

    for element in popList:
        thiosInCrystal.pop(element)

    P.figure()
    n, bins, patches = P.hist(thiosInCrystal, 10, normed=1)
    P.xlabel('dhkl')
    P.ylabel('Occurence')
    P.xlim(0, 7)
    P.savefig('./outputGraphs/'+name+'_smallDist.png')



    print "\n\n----==== DISTRIBUTION RESULTS ====----"
    print "There are", len(thiosInCrystal), "thiophenes with neighbours within 10 angstroms, out of", len(data), "thiophenes considered, representing a portion of", str(float(len(thiosInCrystal))/float(len(data)))+"."


    crystalArray = np.array(thiosInCrystal)
    crystalMean = np.mean(crystalArray)
    crystalStandardDev = np.std(crystalArray)


    print "Only considering the crystalline bit of the distribution, the mean dhkl value =", crystalMean, "and the standard deviation (sqrt of variance) =", str(crystalStandardDev)+". Ideally the mean is as close to 3.76 Angstroms as possible."

    g = crystalStandardDev/float(crystalMean)

    print "Using s^{2} = d^{2} * g^{2}, just using this part of the distribution gives a paracrystallinity of g =", g
    print "---======---"


    print "Now plotting this part of the distribution..."





    return g

#    P.show()



def plotGraph(data, filename, atom1ToHighlight=None, atom2ToHighlight=None):
    thios = []
    alk1 = []
    alk2 = []
    print "Atom 1", atom1ToHighlight
    print "Atom 2", atom2ToHighlight
    for atom in data:
        if atom[1] == 1:
            thios.append([atom[3], atom[4], atom[5]])
        elif atom[1] == 2:
            alk1.append([atom[3], atom[4], atom[5]])
        else:
            alk2.append([atom[3], atom[4], atom[5]])
    fig = P.figure()
    ax = p3.Axes3D(fig)
    thioX, thioY, thioZ = zip(*thios)
    alk1X, alk1Y, alk1Z = zip(*alk1)
    alk2X, alk2Y, alk2Z = zip(*alk2)

    if (atom1ToHighlight == None) or (atom2ToHighlight == None):
        ax.scatter(thioX, thioY, thioZ, s = 20, c = 'r')
        ax.scatter(alk1X, alk1Y, alk1Z, s = 20, c = 'y')
        ax.scatter(alk2X, alk2Y, alk2Z, s = 20, c = 'w')

    else:
        highlightX = [atom1ToHighlight[3], atom2ToHighlight[4][0]]
        highlightY = [atom1ToHighlight[4], atom2ToHighlight[4][1]]
        highlightZ = [atom1ToHighlight[5], atom2ToHighlight[4][2]]
        ax.scatter(thioX, thioY, thioZ, s = 2, c = 'r')
#        ax.scatter(alk1X, alk1Y, alk1Z, s = 2, c = 'y')
#        ax.scatter(alk2X, alk2Y, alk2Z, s = 2, c = 'w')
        #### INCORRECT COORDINATES ####
        ax.scatter(highlightX, highlightY, highlightZ, s = 20, c = 'y')



    ax.set_xlim3d(-70, 70)
    ax.set_ylim3d(-70, 70)
    ax.set_zlim3d(-70, 70)

    P.savefig('./'+filename+'.png')
    


def plotGraph2(data, improper, closestAtom, normalVector):
    thios = []
    alk1 = []
    alk2 = []
    atomPlots = []
    for atom in data:
        if (atom[0] >= improper[0][0]) and (atom[0] <= improper[0][0]+2):
            if atom[1] == 1:
                thios.append([atom[3], atom[4], atom[5]])
            elif atom[1] == 2:
                alk1.append([atom[3], atom[4], atom[5]])
            else:
                alk2.append([atom[3], atom[4], atom[5]])
    fig = P.figure()
    ax = p3.Axes3D(fig)
    thioX, thioY, thioZ = zip(*thios)
    alk1X, alk1Y, alk1Z = zip(*alk1)
    alk2X, alk2Y, alk2Z = zip(*alk2)

    ax.scatter(thioX, thioY, thioZ, s = 20, c = 'r')
    ax.scatter(alk1X, alk1Y, alk1Z, s = 20, c = 'y')
    ax.scatter(alk2X, alk2Y, alk2Z, s = 20, c = 'w')

    piStackingLine = [[improper[0][3]-(30*normalVector[0]), improper[0][4]-(30*normalVector[1]), improper[0][5]-(30*normalVector[2])],
                      [improper[0][3]+(30*normalVector[0]), improper[0][4]+(30*normalVector[1]), improper[0][5]+(30*normalVector[2])]]

    piX, piY, piZ = zip(*piStackingLine)

    highlightX = [improper[0][3], closestAtom[4][0]]
    highlightY = [improper[0][4], closestAtom[4][1]]
    highlightZ = [improper[0][5], closestAtom[4][2]]
    # ax.scatter(thioX, thioY, thioZ, s = 20, c = 'r')
    # ax.scatter(alk1X, alk1Y, alk1Z, s = 20, c = 'y')
    # ax.scatter(alk2X, alk2Y, alk2Z, s = 20, c = 'w')
        #### INCORRECT COORDINATES ####
    ax.scatter(highlightX, highlightY, highlightZ, s = 20, c = 'b')
    ax.plot(piX, piY, piZ, c='k')



    ax.set_xlim3d(-70, 70)
    ax.set_ylim3d(-70, 70)
    ax.set_zlim3d(-70, 70)

    P.savefig('./test.png')
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dat", help="The data file *.dat to run the paracrystallinity measurements on")
    parser.add_argument("-l", "--lammpstrj", help="The corresponding *.lammpstrj file that includes all of the atomic trajectories to run the paracrystallinity measurements on")
    parser.add_argument("-t", "--dumpstep", default=-1, help="The dumpstep within the LAMMSTRJ file desired to obtain the paracrystallinity of. Note that this is not the same as the LAMMPS timestep - it is purely the # of the output within the LAMMPSTRJ file.")
    parser.add_argument("-c", "--cutoff", help="The cut-off in angstroms, beyond which two adjacent thiophenes are said to be in the amorphous regime and before which two adjacent thiophenes are members of the same crystal.")
    parser.add_argument("-p", "--pickle", default=0, help="Is the pickle file containing the crystal pair information already present in the ./outputFiles/ directory? If '1', load the pickle. If '0', then recalculate the crystal pair information.")
    args = parser.parse_args()


    # HARDCODED CHAIN LENGTHS ADDED TO THE MIXES
    # Mw37Additions = [456, 568]
    # Mw46Additions = [246, 540, 544, 532]
    # Mw51Additions = [568, 502, 20, 246, 648, 456, 144, 200]
    # Mw57Additions = [568, 232, 200, 456, 296, 648, 402, 4, 14, 154, 146, 44, 440, 46, 266, 188, 48, 592, 42]

    # CHAIN IDs CORRESPONDING TO THESE LENGTHS:
    Mw37Additions = [2, 1]
    Mw46Additions = [2, 3, 4, 1]
    Mw51Additions = [2, 1, 3, 7, 6, 8, 4, 5]
    Mw57Additions = [19, 17, 10, 11, 9, 13, 7, 1, 44, 5, 20, 14, 2, 12, 6, 18, 3, 15, 4]



#    dumpstepsToBeTreated = np.arange(0, 501, 50)
    # Every dumpstep is at 2,000 lammps timesteps
#    dumpstepsToBeTreated = [0]

    # dumpstepsToBeTreated = np.arange(0, 501, 50)
    paraData = []
    paraDataSmallDist = []
    periodicityOrder = 3 # The maximum number of periodic orders to consider before setting dhkl to an artificial maximum based on the simulation volume
    hkl = [0, 1, 0]

    trajectoryNameRaw = args.lammpstrj
    datNameRaw = args.dat
    dumpstep = int(args.dumpstep)
    crystalCutOff = float(args.cutoff)
    picklePresent = bool(int(args.pickle))

    name = trajectoryNameRaw[:-14]

    if int(name[3]+name[4]) == 37:
        chainLengthAdditions = Mw37Additions
    elif int(name[3]+name[4]) == 46:
        chainLengthAdditions = Mw46Additions
    elif int(name[3]+name[4]) == 51:
        chainLengthAdditions = Mw51Additions
    elif int(name[3]+name[4]) == 57:
        chainLengthAdditions = Mw57Additions
    else:
        chainLengthAdditions = None


    direc = './'
    tolerance = np.pi/6.
    border = 20
    limitThiophenesExamined = False
    blockAdjacent = True
    thiopheneNumber = 100

    print "Tolerance =", tolerance

    trajectoryName = direc+trajectoryNameRaw
    datName = direc+datNameRaw


    trajectoryData, allTimestepNos, simVolData = loadTrajectory(trajectoryName)
    masses, atoms, bonds, angles, dihedrals, impropers = loadDat(datName)

    # TO TRY AND GET THE RIGHT MOLECULES FOR THE LENGTH
    # chainAtoms = {}
    # for atom in atoms:
    #     if atom[1] not in chainAtoms:
    #         chainAtoms[atom[1]] = [atom]
    #         continue
    #     chainAtoms[atom[1]].append(atom)
    # for chain in chainAtoms:
    #     print chain, len(chainAtoms[chain])/3.
    # raise SystemError('STOP')

    timestepNos = []

    if dumpstep >= len(trajectoryData):
        # Just run the final dumpstep
        dumpstep = -1


    # periodicTraj = generatePeriodicity(trajectoryData[dumpstep]. simVolData[dumpstep])



    print "Treating dumpstep", str(dumpstep)+"..."

    outputDirectory = './outputFiles/'
    crystalPairsPickle = name+'_crystalPairs.pickle'

    if picklePresent == True:
        print "Pickle found, loading crystalpairs..."
        pickleFile = open(outputDirectory+crystalPairsPickle, 'r')
        crystalPairs = pickle.load(pickleFile)
    else:
        print "Pickle not found, calculating crystalpairs..."
        paraG, dhkl, crystalPairs = calculateCrystals(atoms, bonds, angles, dihedrals, impropers, trajectoryData[dumpstep], simVolData[dumpstep], hkl, periodicityOrder, tolerance, crystalCutOff, border, limitThiophenesExamined, thiopheneNumber, blockAdjacent)
        print "Dumping crystalpairs as pickle..."
        pickleFile = open(outputDirectory+crystalPairsPickle, 'w+')
        pickle.dump(crystalPairs, pickleFile)

    print "There are", len(crystalPairs), "crystal pairs."


    bondDict = getBondDict(bonds)
    crystals, crystalsAssorted = treatCrystalPairs(crystalPairs, bondDict)



    crystalsPerChain = getCrystalsPerChain(atoms, crystals)
    chainLengths = getChainLengths(atoms)



    # for wobbey in crystalsPerChain:
    #     print "Chain =", wobbey, "has crystals numbering", len(crystalsPerChain[wobbey])


    # Examine chain 7

    # print crystalsPerChain[7]
    # for crystalID in crystalsPerChain[7]:
    #     print "CrystalID = ", crystalID, "Crystal atoms =", crystals[crystalID]


    # raise SystemError('STOP')

    # makeLammpsFileChainTest(name, 7, crystalsPerChain[7], crystals)

    makeCSVs(name, crystalsPerChain, crystals, chainLengths)

    biggestCrystal = getBiggestCrystal(crystals)
    makeNewLammpsFiles(name, crystals, biggestCrystal, chainLengthAdditions)

    print "New Lammps file created."









    # crystalPairs = [[7, 319], [10, 6136], [10, 319], [25, 5626], [25, 3586], [31, 3598], [34, 3598], [43, 6124], [43, 223], [49, 151], [49, 6121], [49, 6118], [52, 148], [73, 5863], [79, 5860], [82, 5857], [97, 5848], [100, 5848], [103, 5845], [124, 5953], [133, 241], [136, 238], [139, 235], [151, 214], [52, 151], [163, 283], [166, 283], [184, 5863], [190, 5872], [193, 262], [193, 265], [208, 418], [208, 415], [46, 220], [46, 223], [229, 400], [232, 397], [232, 400], [235, 397], [241, 385], [244, 382], [244, 385], [130, 244], [256, 5953], [259, 607], [259, 604], [196, 262], [190, 265], [265, 487], [283, 355], [283, 352], [160, 289], [10, 316], [328, 11368], [328, 532], [331, 535], [334, 11371], [334, 538], [337, 11374], [337, 11377], [340, 11377], [343, 11380], [349, 433], [352, 430], [355, 430], [355, 427], [211, 364], [370, 412], [376, 877], [382, 880], [382, 883], [385, 883], [238, 391], [235, 394], [397, 3619], [400, 3616], [403, 3613], [370, 415], [205, 415], [418, 496], [361, 421], [358, 424], [349, 436], [439, 4093], [442, 4096], [451, 11278], [457, 11272], [457, 6532], [460, 11269], [460, 6529], [460, 11266], [100, 469], [472, 586], [472, 583], [475, 721], [478, 724], [481, 727], [481, 730], [265, 484], [490, 613], [493, 616], [376, 499], [370, 502], [505, 661], [508, 658], [511, 655], [520, 643], [328, 535], [553, 1048], [556, 1045], [556, 4093], [559, 4093], [559, 1042], [559, 4096], [562, 1039], [568, 4102], [469, 583], [586, 718], [487, 610], [616, 742], [427, 625], [424, 625], [628, 1060], [634, 1159], [640, 766], [523, 640], [526, 640], [523, 643], [517, 646], [661, 868], [664, 871], [667, 8977], [667, 874], [670, 877], [673, 883], [673, 880], [676, 886], [679, 892], [682, 892], [682, 10465], [694, 7831], [697, 7828], [703, 8737], [709, 8731], [589, 715], [718, 8725], [718, 8722], [721, 8722], [613, 739], [616, 739], [742, 916], [751, 868], [760, 1168], [766, 1249], [643, 766], [769, 1246], [775, 826], [535, 775], [799, 8023], [802, 8020], [808, 8026], [811, 8029], [826, 7975], [835, 7873], [841, 1240], [841, 9844], [841, 9841], [844, 9838], [847, 1177], [847, 1174], [847, 9835], [847, 9838], [751, 865], [667, 877], [379, 880], [913, 967], [916, 970], [916, 967], [745, 916], [919, 973], [925, 979], [928, 1072], [928, 982], [931, 1075], [931, 985], [934, 988], [934, 1075], [934, 1078], [949, 1015], [952, 1018], [967, 1444], [970, 1447], [919, 970], [922, 973], [973, 1450], [922, 976], [925, 976], [979, 1339], [988, 1330], [991, 1327], [937, 991], [994, 1324], [997, 1321], [1009, 1105], [1012, 8824], [952, 1015], [1015, 8824], [568, 1033], [565, 1036], [1048, 1144], [1051, 1150], [1054, 1153], [928, 1069], [1033, 1069], [925, 1069], [931, 1072], [937, 1078], [1081, 1117], [1084, 1120], [1084, 1123], [1099, 8812], [1102, 8815], [1105, 8818], [1006, 1105], [1111, 8713], [1114, 8713], [1129, 4075], [1132, 4078], [1132, 4081], [1138, 1273], [1150, 1261], [1057, 1156], [1171, 1207], [850, 1174], [1216, 9580], [1219, 9577], [1222, 9724], [1225, 9574], [1225, 9721], [1225, 9724], [1228, 9571], [1183, 1231], [844, 1240], [766, 1246], [1156, 1255], [1153, 1261], [1261, 1351], [1147, 1264], [1267, 7996], [1270, 7999], [1273, 8002], [1285, 1375], [1288, 1378], [1291, 8479], [1294, 8482], [1297, 8485], [1309, 1393], [1312, 1393], [1000, 1321], [1324, 9361], [985, 1333], [1333, 1615], [982, 1336], [1348, 1462], [1354, 7993], [1357, 7996], [1357, 1471], [1360, 7999], [1360, 7996], [1360, 1543], [1363, 8002], [1369, 1552], [1372, 1555], [1375, 1558], [1378, 1561], [1381, 1564], [1291, 1381], [1294, 1381], [1384, 1567], [1414, 1633], [1414, 1636], [1417, 1633], [1420, 1630], [949, 1423], [955, 1429], [961, 1435], [1438, 9313], [964, 1444], [1456, 1525], [1456, 1528], [1459, 1528], [1459, 1531], [1462, 1531], [1465, 1534], [1351, 1465], [1354, 1468], [1483, 1687], [1486, 1684], [1489, 1597], [1501, 1609], [1504, 9286], [1504, 1609], [1504, 1612], [1516, 9343], [1519, 9343], [1522, 9340], [1525, 9337], [1525, 9334], [1534, 9433], [1468, 1537], [1363, 1546], [1366, 1549], [1552, 1603], [1378, 1558], [1561, 9370], [1570, 1654], [1573, 1657], [1576, 1660], [1579, 7507], [1588, 1669], [1492, 1597], [1492, 1600], [1549, 1603], [1507, 1612], [1510, 1615], [1507, 1615], [1333, 1618], [1621, 9352], [1645, 3679], [1645, 3676], [1648, 3673], [1651, 1708], [1657, 5689], [1585, 1666], [1672, 10918], [1675, 10921], [1675, 10918], [1678, 10924], [1483, 1684], [1648, 1708], [1711, 3673], [1726, 10027], [1732, 4795], [1735, 4798], [1741, 1936], [1753, 1846], [1753, 1849], [1756, 1849], [1759, 1990], [1774, 2005], [1777, 2008], [1780, 6826], [1780, 2008], [1792, 2038], [1792, 2035], [1795, 2038], [1810, 4852], [1810, 7345], [1819, 7204], [1819, 7207], [1822, 7201], [1825, 7201], [1828, 2071], [1834, 7375], [1837, 7375], [1846, 1981], [1852, 1924], [1858, 3730], [1858, 9145], [1861, 3736], [1861, 9148], [1867, 9154], [1876, 2602], [1879, 2599], [1882, 2011], [1882, 2014], [1882, 2596], [1885, 2011], [1888, 2158], [1891, 2161], [1900, 2254], [1903, 2257], [1849, 1921], [1855, 1924], [1933, 1972], [1738, 1936], [1945, 4576], [1951, 2413], [1954, 2413], [1966, 2302], [1975, 2107], [1978, 2104], [1843, 1984], [1918, 1987], [1915, 1990], [1765, 1993], [1762, 1993], [1768, 1999], [1771, 2002], [1777, 2005], [1795, 2041], [2062, 2344], [2062, 2341], [1828, 2068], [2071, 2356], [2074, 2356], [2077, 7384], [2083, 4549], [2083, 4546], [2086, 4552], [2101, 2314], [2101, 2317], [1975, 2104], [2113, 2293], [2116, 2290], [2116, 3718], [2119, 3715], [2122, 11389], [2128, 3694], [2134, 2617], [2137, 2614], [2140, 2272], [2143, 2269], [2143, 2458], [2143, 2272], [2146, 2458], [2146, 2266], [2146, 2461], [2149, 2461], [2152, 2464], [2164, 2581], [2182, 8899], [2185, 8896], [2185, 8893], [2191, 4867], [2194, 4867], [2194, 2230], [2209, 7315], [2209, 7318], [2242, 2482], [2254, 2473], [2257, 2470], [1903, 2260], [1906, 2260], [2137, 2275], [2278, 3733], [2278, 3730], [2284, 2443], [1966, 2299], [2299, 2425], [2302, 2422], [2098, 2314], [2320, 2542], [2323, 2545], [2326, 2548], [2326, 2545], [2329, 2551], [2338, 2488], [2059, 2341], [2341, 2488], [2341, 2491], [2344, 2491], [2068, 2350], [2071, 2353], [2068, 2353], [2359, 2503], [2377, 2530], [2377, 4558], [2380, 4561], [2383, 4564], [2389, 4324], [2395, 4609], [2401, 5230], [2404, 5227], [2416, 5242], [2419, 5245], [2419, 5248], [2422, 5248], [2290, 2437], [2449, 2785], [2461, 2656], [2461, 2659], [2155, 2467], [2254, 2476], [2251, 2476], [2479, 2566], [2251, 2479], [2479, 2563], [2356, 2500], [2500, 2995], [2362, 2506], [2518, 4306], [2521, 4309], [2377, 2527], [2380, 2530], [2542, 2830], [2545, 2833], [2332, 2554], [2560, 2851], [2569, 2941], [2584, 7855], [2158, 2587], [2593, 5530], [2593, 5533], [2134, 2614], [2635, 4243], [2608, 2644], [2647, 11332], [2608, 2647], [2683, 2770], [2686, 2767], [2689, 2764], [2527, 2701], [2704, 11077], [2707, 11083], [2713, 11089], [2719, 11095], [2734, 6166], [2743, 2989], [2746, 2986], [2749, 2983], [2761, 2887], [2764, 2884], [2686, 2764], [2770, 2818], [2770, 2821], [2449, 2782], [2800, 9895], [2803, 9898], [2806, 8935], [2809, 8938], [2815, 2878], [2818, 2875], [2680, 2824], [2680, 2827], [2842, 3001], [2845, 3001], [2848, 3004], [2560, 2854], [2557, 2854], [2863, 2962], [2671, 2863], [2866, 2962], [2869, 2965], [2875, 3028], [2818, 2878], [2878, 3031], [2893, 10945], [2896, 10948], [2899, 10951], [2905, 2983], [2908, 2986], [2908, 3067], [2911, 3067], [2917, 2998], [2920, 6091], [2920, 3001], [2950, 7876], [2953, 7876], [2953, 7879], [2956, 7879], [2968, 3022], [2911, 2989], [2914, 2992], [2914, 2998], [2917, 3001], [3016, 7888], [3016, 7891], [3025, 7903], [3028, 7906], [3031, 7909], [2875, 3031], [2878, 3034], [2881, 3034], [3046, 10945], [3052, 10951], [2908, 3064], [2905, 3064], [3067, 6172], [3076, 9619], [3085, 9880], [3097, 11440], [3103, 11437], [3103, 11434], [3115, 10150], [3118, 10153], [3118, 3250], [3121, 10156], [3124, 3256], [3124, 10159], [3127, 10282], [3142, 3271], [3142, 6589], [3145, 3274], [3151, 3280], [3154, 10138], [3157, 3286], [3157, 10135], [3160, 3289], [3166, 4207], [3166, 3295], [3172, 3304], [3175, 3307], [3178, 3310], [3181, 3313], [3181, 3316], [3190, 9829], [3205, 9844], [3208, 9841], [3214, 5419], [3217, 5419], [3232, 3328], [3235, 3331], [3235, 3334], [3238, 3334], [3244, 3340], [3121, 3253], [3283, 3445], [3160, 3292], [3163, 3292], [3238, 3295], [3169, 3301], [3100, 3301], [3097, 3304], [3178, 3307], [3094, 3307], [3181, 3310], [3088, 3310], [3184, 3313], [3316, 9736], [3319, 9739], [3241, 3337], [3247, 3340], [3361, 3589], [3361, 10210], [3370, 10225], [3373, 10228], [3385, 10237], [3388, 10240], [3391, 10243], [3391, 6139], [3394, 10246], [3397, 10249], [3397, 3574], [3400, 3571], [3400, 10252], [3403, 3568], [3427, 3535], [3430, 3532], [3433, 3529], [3340, 3433], [3436, 3526], [3337, 3436], [3454, 3853], [3457, 3856], [3460, 3856], [3463, 3634], [3463, 3631], [3481, 7129], [3493, 7150], [3496, 7147], [3496, 7144], [3496, 6784], [3499, 6790], [3502, 6793], [3514, 9751], [3526, 4036], [3529, 3625], [3532, 3625], [3547, 3658], [3565, 6112], [3577, 6130], [3577, 6133], [3358, 3589], [25, 3589], [3601, 5425], [3604, 5422], [400, 3613], [3526, 3619], [394, 3619], [394, 3622], [3466, 3637], [3667, 6118], [3685, 3754], [3688, 3751], [2128, 3697], [2125, 3700], [3703, 11389], [2113, 3721], [1855, 3730], [1858, 3733], [3748, 9133], [3688, 3748], [3685, 3751], [3769, 9703], [3772, 9706], [3775, 9709], [3787, 9217], [3805, 5284], [3811, 3871], [3811, 3868], [3814, 3874], [3817, 3877], [3829, 4834], [3829, 3943], [3829, 3940], [3832, 4837], [3835, 3937], [3835, 4840], [3838, 3934], [3838, 4843], [3838, 3931], [3841, 4846], [3451, 3850], [3853, 5686], [3859, 3991], [3862, 3994], [3868, 5200], [3871, 5203], [3814, 3871], [3874, 5206], [3877, 5209], [3880, 5212], [3883, 5215], [3886, 5218], [3886, 5215], [3895, 4837], [3895, 4834], [3898, 4840], [3901, 4006], [3901, 4003], [3901, 4843], [3901, 4840], [3904, 4843], [3904, 4846], [3907, 4846], [3910, 4018], [3835, 3934], [3823, 3946], [3820, 3949], [3817, 3952], [3814, 3955], [3955, 6247], [3955, 6244], [3958, 6250], [3472, 3979], [3472, 3982], [3862, 3991], [3994, 5191], [3997, 5188], [4000, 5185], [3898, 4003], [4006, 5674], [4009, 5674], [4009, 5671], [4018, 7630], [4027, 7645], [4030, 7645], [4042, 4117], [4048, 4081], [4048, 4084], [4045, 4084], [556, 4090], [562, 4093], [4108, 5326], [4111, 7597], [4111, 7600], [4114, 5323], [4045, 4120], [4123, 4216], [4123, 4213], [4126, 4165], [4129, 4168], [4135, 11212], [4138, 11206], [4141, 10732], [4141, 11203], [4156, 11233], [4156, 11236], [4159, 11230], [4159, 11233], [4162, 11227], [4165, 11224], [4168, 11221], [4126, 4168], [4171, 8110], [4189, 9784], [4192, 9880], [4192, 9781], [4192, 9784], [4195, 9781], [4195, 9883], [4198, 9775], [4198, 9778], [3169, 4204], [2632, 4246], [4261, 4483], [4264, 11038], [4267, 11044], [4270, 11308], [4273, 11305], [4279, 11128], [4285, 11125], [4285, 4375], [4285, 11122], [4291, 11116], [4291, 11113], [4294, 4384], [4294, 4387], [4297, 4387], [4300, 10981], [4300, 10978], [4303, 10981], [4312, 10990], [2524, 4312], [4315, 10993], [4318, 10996], [4327, 4564], [4330, 4564], [4330, 4603], [4333, 4600], [4336, 4597], [4339, 4591], [4375, 10996], [4288, 4378], [4291, 4381], [301, 4396], [4402, 4648], [4411, 4639], [4414, 4639], [4417, 11134], [4417, 4633], [4423, 4486], [4426, 4483], [4429, 4624], [4429, 11014], [4438, 5065], [4453, 4834], [4468, 7534], [4261, 4486], [4264, 4486], [4492, 11143], [4495, 11146], [4498, 11149], [4504, 8671], [4510, 11338], [4510, 8845], [4516, 11344], [4522, 4678], [4528, 4675], [4537, 7393], [4540, 4654], [2080, 4546], [1945, 4573], [1942, 4576], [4582, 4783], [4342, 4591], [4591, 7267], [4339, 4594], [4609, 4732], [4609, 4735], [2395, 4612], [4426, 4627], [4420, 4633], [4633, 4702], [4636, 4699], [4669, 7297], [4531, 4672], [4678, 7246], [4525, 4678], [4522, 4681], [4639, 4696], [4696, 4750], [4696, 4753], [4639, 4699], [4630, 4705], [4627, 4708], [4708, 7807], [4708, 7810], [4711, 7807], [4714, 7804], [4615, 4726], [4612, 4729], [4606, 4735], [4699, 4750], [4756, 7822], [4759, 7255], [4768, 7102], [4774, 7264], [4774, 7783], [4774, 7780], [4777, 7783], [4777, 7786], [4786, 7735], [4789, 7738], [1735, 4795], [4795, 5743], [1738, 4798], [4819, 7411], [3892, 4834], [3835, 4837], [1807, 4855], [4873, 8890], [4873, 8893], [4882, 8968], [4885, 6733], [4888, 6730], [4891, 6727], [4891, 4981], [4891, 4978], [4894, 6724], [4894, 5110], [4897, 5107], [4897, 6721], [4897, 5104], [4906, 6709], [4906, 6712], [4909, 6706], [4909, 6709], [4909, 5206], [4912, 6703], [4912, 6706], [4921, 5191], [4921, 5188], [4927, 6781], [4930, 6778], [4933, 5077], [4933, 6775], [4933, 6772], [4936, 5077], [4936, 6772], [4945, 10567], [4948, 5059], [4948, 10570], [4948, 10573], [4951, 5056], [4951, 10573], [4954, 10579], [4969, 6748], [4975, 5155], [4978, 5152], [4993, 8194], [4999, 8188], [5005, 5575], [5014, 10108], [5014, 10111], [5017, 5140], [5017, 10111], [5029, 5149], [5032, 5152], [5035, 5155], [5035, 8854], [4957, 5050], [4954, 5053], [5059, 10672], [4945, 5059], [5062, 10669], [4441, 5065], [4438, 5068], [4936, 5080], [5083, 5179], [5086, 5176], [5092, 5170], [5098, 5167], [5098, 5164], [5101, 5161], [4894, 5107], [5113, 5272], [5113, 5275], [5119, 8197], [5122, 8194], [5128, 8188], [5014, 5140], [5020, 5143], [5023, 5146], [4981, 5152], [4978, 5155], [4975, 5158], [5038, 5158], [5101, 5164], [5095, 5167], [5173, 5218], [5089, 5173], [5173, 5221], [4003, 5185], [3865, 5197], [4912, 5206], [5227, 10681], [5230, 10681], [5230, 10684], [5233, 10684], [5236, 10687], [5026, 5257], [5026, 5260], [5260, 10069], [5260, 10072], [3808, 5284], [5296, 6229], [5314, 6490], [5317, 6490], [3931, 5350], [3934, 5353], [3931, 5353], [5371, 7420], [5374, 9517], [5392, 5533], [5404, 5512], [5404, 5515], [5407, 5509], [3211, 5413], [5416, 5620], [5434, 11437], [5449, 8566], [5452, 5596], [5452, 8569], [5470, 5569], [5473, 5572], [5476, 5578], [2029, 5476], [5479, 5581], [5482, 5584], [5485, 5587], [5485, 5584], [5485, 5644], [5488, 5590], [5491, 5593], [5494, 5596], [5494, 5635], [5503, 5629], [5515, 11356], [5527, 11347], [5527, 11344], [5542, 7855], [5545, 7855], [5560, 6475], [5563, 6478], [5566, 6481], [5569, 6484], [5008, 5572], [5008, 5575], [5581, 8581], [5584, 8581], [5452, 5593], [5497, 5599], [5602, 8560], [28, 5623], [5506, 5623], [5503, 5626], [28, 5626], [5491, 5638], [5488, 5641], [5647, 10480], [5650, 10477], [5650, 10480], [5662, 10468], [5665, 10462], [5665, 10465], [5680, 7636], [5683, 7639], [5692, 9523], [5692, 9520], [1717, 5692], [5383, 5695], [5380, 5695], [5701, 9529], [5710, 5746], [5713, 5749], [5716, 5893], [5719, 5893], [5719, 5890], [5722, 5890], [5722, 5887], [5734, 7735], [5713, 5746], [5752, 7750], [5755, 7753], [5764, 6955], [5767, 7768], [5767, 6952], [5785, 5935], [5785, 5989], [5788, 5992], [5794, 5905], [5797, 9550], [5800, 6955], [5803, 6952], [5806, 6949], [5809, 6946], [5809, 6949], [5812, 6196], [5824, 6598], [5824, 6601], [5842, 5881], [5854, 10129], [5857, 10132], [76, 5863], [193, 5875], [5878, 5944], [5881, 5941], [5881, 5938], [5848, 5884], [4807, 5896], [5794, 5902], [5911, 9556], [5914, 9559], [5923, 5980], [5785, 5932], [5782, 5935], [5884, 5938], [5878, 5941], [5833, 5944], [121, 5950], [5917, 5983], [5785, 5992], [6001, 9472], [4807, 6001], [4810, 6001], [6013, 9382], [1684, 6022], [6031, 9451], [6031, 10957], [6031, 10960], [6034, 9457], [6034, 9454], [6049, 10978], [6058, 6154], [6070, 6184], [6070, 6181], [5974, 6070], [5977, 6070], [6073, 6184], [6079, 6298], [3565, 6109], [6112, 6154], [46, 6124], [3574, 6133], [6061, 6157], [2737, 6163], [6172, 9406], [3064, 6175], [5812, 6199], [6208, 6268], [3637, 6211], [3634, 6214], [3634, 6217], [6226, 10312], [5296, 6232], [3952, 6241], [6247, 8218], [6250, 8215], [3547, 6280], [6280, 10177], [3550, 6283], [3553, 6283], [6286, 10183], [3553, 6286], [3556, 6286], [6082, 6295], [6307, 9229], [6307, 9226], [6310, 9226], [6316, 9019], [6328, 8143], [6340, 8359], [6343, 8356], [6346, 8353], [6370, 6427], [6379, 6442], [6385, 10183], [6388, 10372], [6391, 10375], [6403, 10333], [6406, 10387], [6412, 10396], [6415, 8425], [6436, 9241], [6439, 9241], [6439, 9244], [6442, 9244], [6442, 9247], [6493, 6541], [6496, 6586], [6502, 6580], [6511, 7516], [6514, 7519], [6514, 7516], [6517, 7522], [6520, 7522], [6496, 6538], [4231, 6541], [6550, 11251], [6550, 11254], [6559, 10966], [6565, 9466], [6568, 9469], [5368, 6574], [5365, 6577], [3793, 6607], [6625, 6970], [6625, 6973], [6631, 7360], [6634, 7357], [6634, 6814], [6637, 7354], [6637, 6814], [6643, 7348], [6646, 6880], [6649, 6880], [6652, 7069], [6652, 7072], [6670, 7279], [6670, 7282], [6673, 7759], [6676, 7762], [6679, 7093], [6679, 7765], [6685, 7090], [6685, 7087], [6688, 7087], [6700, 7174], [4891, 6724], [6724, 7345], [4882, 6736], [4876, 6739], [4873, 6742], [6751, 7327], [6751, 7330], [6754, 7324], [6775, 7162], [6778, 7159], [3493, 6784], [3496, 6787], [3502, 6790], [6793, 7639], [6658, 6802], [6814, 6982], [6817, 6985], [6823, 6919], [6826, 6919], [6850, 8170], [6853, 7003], [6853, 8170], [6859, 6994], [6859, 6997], [6862, 6991], [6862, 6994], [6865, 6988], [6871, 8200], [6871, 8203], [3967, 6889], [6898, 7036], [6901, 7033], [6904, 7030], [6907, 7027], [6910, 7024], [6913, 7021], [6862, 6913], [6916, 7018], [6823, 6922], [6928, 6979], [6928, 8131], [6934, 7042], [6937, 7042], [6610, 6946], [5764, 6952], [5761, 6958], [6967, 7054], [6628, 6973], [6931, 6976], [6865, 6991], [1792, 6994], [6856, 7000], [6850, 7006], [6847, 7006], [7021, 8158], [7024, 8155], [7027, 8152], [7030, 8149], [7033, 8146], [6898, 7033], [7036, 8143], [7036, 8140], [7042, 8137], [6934, 7045], [6931, 7045], [6970, 7051], [7072, 7210], [7081, 7180], [7081, 7177], [7087, 7291], [7087, 7294], [7090, 7291], [6682, 7093], [7108, 7258], [4351, 7108], [7114, 7252], [7117, 7252], [3481, 7132], [3484, 7132], [3499, 7144], [7153, 10480], [6781, 7159], [6778, 7162], [7168, 7228], [6697, 7177], [7084, 7183], [7189, 7285], [1825, 7198], [7069, 7210], [7168, 7225], [7165, 7228], [7228, 7309], [4681, 7246], [4594, 7267], [6673, 7279], [7192, 7285], [7090, 7288], [4669, 7300], [4672, 7300], [7312, 10450], [7315, 10447], [7318, 10444], [6754, 7327], [6748, 7330], [4864, 7336], [6721, 7342], [4849, 7354], [1837, 7372], [4537, 7390], [7396, 10243], [7396, 10246], [7405, 10447], [7417, 9511], [5368, 7417], [7420, 9511], [7435, 10405], [7450, 8878], [7468, 8911], [7471, 8914], [7480, 9799], [7483, 9799], [7489, 8461], [7489, 8464], [7495, 8542], [7495, 8545], [7501, 8539], [7501, 8536], [7504, 8536], [1582, 7507], [6517, 7519], [4465, 7534], [4468, 7537], [7537, 8701], [7549, 7696], [7552, 7699], [7561, 7624], [7564, 7621], [7576, 7675], [7582, 8710], [7582, 8707], [7585, 8713], [5335, 7585], [7588, 8713], [7588, 8716], [4108, 7594], [4108, 7597], [7609, 7666], [4015, 7627], [6790, 7636], [6796, 7642], [7651, 9754], [7657, 9763], [7660, 9766], [7684, 8809], [7546, 7693], [7543, 7693], [7546, 7696], [7699, 8587], [7717, 7801], [7720, 7798], [7723, 7798], [5737, 7735], [5743, 7741], [5746, 7744], [5749, 7747], [5755, 7750], [5758, 7753], [5764, 7765], [4759, 7825], [697, 7831], [691, 7834], [7837, 8749], [7840, 8752], [2581, 7855], [2584, 7858], [2950, 7873], [3016, 7894], [3019, 7897], [3022, 7897], [7909, 9565], [7918, 9907], [7921, 9907], [7945, 11302], [7948, 11305], [7948, 8056], [7948, 11302], [7951, 8053], [7951, 11308], [7954, 8050], [7954, 11311], [7957, 8047], [7957, 11314], [823, 7978], [820, 7981], [1357, 7993], [1366, 8002], [1276, 8005], [8011, 8092], [799, 8020], [8023, 8068], [802, 8023], [814, 8032], [8035, 9895], [8026, 8071], [7918, 8074], [7912, 8077], [7915, 8077], [8014, 8095], [8014, 8098], [4174, 8107], [8122, 9007], [8125, 9007], [6322, 8140], [7039, 8140], [6325, 8140], [6331, 8146], [8152, 8251], [8152, 8254], [8155, 8302], [8167, 8290], [4999, 8185], [5002, 8185], [4996, 8188], [4996, 8191], [4993, 8191], [6250, 8212], [8227, 10015], [8233, 10012], [8239, 8524], [8251, 8374], [8149, 8254], [8254, 8374], [8263, 9070], [6322, 8263], [8293, 8383], [8293, 8380], [8314, 10015], [8320, 10021], [5278, 8323], [7417, 8341], [8344, 9511], [8359, 8419], [3763, 8362], [8254, 8371], [8377, 9604], [8296, 8380], [8392, 8494], [8395, 8488], [8401, 8482], [8410, 9616], [8413, 8536], [8413, 8533], [8362, 8419], [3763, 8422], [3766, 8425], [8455, 9808], [7450, 8455], [8458, 8548], [8458, 9805], [8467, 11467], [8404, 8476], [8488, 9772], [8392, 8491], [8389, 8494], [8497, 9652], [8389, 8497], [8500, 9655], [8509, 9994], [8512, 9991], [8236, 8524], [8410, 8539], [7498, 8542], [7492, 8545], [8563, 8620], [5452, 8566], [8566, 8617], [5455, 8569], [8602, 11026], [8605, 11407], [8569, 8614], [8614, 9931], [8569, 8617], [8617, 9928], [8566, 8620], [8629, 11446], [8629, 11449], [8632, 11449], [8638, 11455], [8647, 9916], [8650, 9913], [8659, 9961], [8602, 8665], [7540, 8698], [1111, 8716], [715, 8725], [712, 8728], [715, 8728], [8731, 8830], [8731, 8833], [8734, 8833], [706, 8734], [8740, 8839], [700, 8740], [8740, 8842], [8743, 8842], [7834, 8746], [8746, 8845], [8746, 8848], [8749, 8848], [4237, 8770], [4240, 8770], [4243, 8776], [4249, 8782], [8794, 9976], [8797, 9979], [8800, 9982], [1096, 8809], [7687, 8809], [7684, 8812], [7681, 8812], [1108, 8818], [8728, 8827], [1015, 8827], [8737, 8836], [8734, 8836], [8737, 8839], [8674, 8839], [8671, 8842], [5032, 8857], [8869, 10600], [8872, 10600], [7453, 8878], [2182, 8896], [6736, 8899], [7465, 8908], [2803, 8935], [664, 8980], [8989, 9322], [8989, 9319], [8278, 8992], [8992, 9322], [8275, 8995], [9001, 9166], [9022, 9229], [9028, 9235], [9028, 9232], [9031, 9235], [9031, 9238], [9049, 9250], [9049, 9247], [9064, 9589], [8263, 9067], [8989, 9082], [9097, 9589], [9103, 9334], [9109, 9175], [9109, 9340], [9112, 9169], [9115, 9166], [9121, 9304], [9124, 9298], [9133, 9196], [8998, 9166], [9109, 9172], [9013, 9178], [9010, 9181], [3784, 9184], [3781, 9187], [7066, 9193], [9202, 9535], [9205, 9289], [9208, 9286], [9214, 9490], [9217, 9490], [9220, 9493], [9025, 9232], [9046, 9247], [9256, 9508], [9259, 9421], [9271, 9496], [9271, 9433], [9274, 9439], [9298, 9349], [9121, 9301], [8992, 9319], [9106, 9337], [9109, 9337], [1519, 9340], [1516, 9346], [9301, 9349], [1618, 9355], [9385, 9481], [9388, 9484], [9394, 9448], [9397, 9451], [9397, 9454], [9256, 9418], [9274, 9436], [6025, 9445], [6028, 9451], [9394, 9454], [6037, 9457], [5998, 9469], [6565, 9469], [9214, 9487], [9274, 9493], [5374, 9520], [5695, 9523], [5695, 9526], [5908, 9553], [1213, 9583], [9067, 9592], [9598, 9688], [9601, 9685], [8374, 9601], [3079, 9619], [9622, 9700], [9625, 9697], [9625, 9700], [9628, 9730], [9631, 9730], [9631, 9733], [9634, 9688], [9634, 9691], [9637, 9736], [9637, 9739], [9640, 9739], [9640, 9742], [8497, 9649], [8494, 9649], [8500, 9652], [9652, 9760], [8503, 9655], [9598, 9685], [9631, 9691], [3769, 9706], [9718, 9817], [9721, 9820], [9721, 9823], [9634, 9733], [9634, 9736], [3322, 9742], [9652, 9757], [4195, 9778], [4066, 9778], [4063, 9781], [9799, 11464], [9817, 9868], [9715, 9817], [9718, 9820], [9820, 9871], [3190, 9826], [844, 9841], [3208, 9844], [3205, 9847], [9820, 9868], [8032, 9895], [2800, 9898], [2806, 9901], [8071, 9901], [7918, 9904], [7921, 9910], [8614, 9928], [9934, 11413], [9937, 11410], [9940, 11407], [9943, 11041], [9943, 11038], [9967, 10054], [9970, 10051], [1402, 9982], [8509, 9991], [8179, 9994], [8179, 9997], [8236, 10009], [8230, 10015], [8227, 10018], [8317, 10018], [1726, 10030], [1960, 10033], [10045, 10120], [10048, 10120], [10048, 10123], [10054, 10114], [9970, 10054], [10057, 10111], [5257, 10066], [5263, 10072], [10081, 10123], [5131, 10081], [5134, 10081], [10084, 10126], [3118, 10150], [3124, 10156], [6283, 10180], [10183, 10375], [10186, 10372], [10189, 10372], [3412, 10192], [10201, 10264], [10204, 10264], [10204, 10267], [3355, 10204], [3358, 10210], [3376, 10231], [3379, 10234], [3382, 10234], [10258, 10354], [10261, 10354], [10261, 10351], [10264, 10351], [10156, 10276], [10276, 10642], [3130, 10285], [10285, 10633], [3136, 10291], [3271, 10312], [10324, 10390], [6403, 10330], [10339, 10423], [10342, 10426], [10351, 10435], [10258, 10357], [2227, 10363], [6391, 10372], [10327, 10390], [6412, 10393], [6415, 10396], [7435, 10402], [10411, 10624], [7324, 10438], [7408, 10444], [5659, 10471], [5656, 10474], [7153, 10477], [5638, 10522], [10534, 11350], [10546, 11362], [10549, 11365], [10510, 10558], [4942, 10567], [4957, 10579], [4957, 10582], [7486, 10612], [10618, 10837], [10624, 10897], [10627, 10846], [10630, 10846], [10282, 10636], [10273, 10642], [10645, 10729], [10273, 10645], [8860, 10696], [10702, 10795], [10702, 10798], [10597, 10708], [10723, 10855], [10642, 10726], [4144, 10732], [10750, 10849], [10750, 10846], [10753, 10849], [10759, 10858], [10762, 10861], [10765, 11197], [10771, 11062], [10771, 11065], [10774, 11065], [5050, 10780], [5047, 10783], [10789, 10873], [10792, 10876], [10792, 10873], [10699, 10795], [10795, 10876], [10837, 10891], [10621, 10840], [10624, 10840], [10627, 10843], [10762, 10858], [10723, 10858], [10873, 11077], [10888, 10936], [10840, 10891], [10891, 10933], [10297, 10903], [1672, 10915], [3049, 10945], [3052, 10948], [3049, 10948], [10924, 10966], [6556, 10966], [8602, 11029], [4261, 11035], [11011, 11053], [11011, 11056], [11098, 11248], [11116, 11179], [11119, 11176], [4417, 11137], [4420, 11137], [4489, 11140], [11143, 11314], [11146, 11317], [11146, 11320], [11149, 11323], [11152, 11326], [11173, 11377], [11173, 11374], [11176, 11380], [11119, 11179], [11188, 11290], [7936, 11221], [7939, 11224], [11233, 11290], [4153, 11236], [11236, 11287], [4234, 11245], [454, 11272], [451, 11275], [454, 11275], [11233, 11287], [11185, 11290], [11191, 11293], [7945, 11299], [4276, 11302], [4270, 11305], [11149, 11326], [11152, 11329], [10537, 11353], [5518, 11356], [331, 11368], [331, 11371], [11170, 11374], [334, 11374], [11176, 11377], [9943, 11404], [8608, 11410], [9940, 11410], [8611, 11413], [11212, 11422], [11212, 11425], [11209, 11425], [3106, 11434], [3100, 11437], [3100, 11440], [8635, 11452], [9790, 11467]] # 2): ('Mw_37_PDI_2.3_con_40_Equil_623K.dat', 'Mw_37_PDI_2.3_con_40_Equil_623K.dat.lammpstrj') Crystal Pairs

