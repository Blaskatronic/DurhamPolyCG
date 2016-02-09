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
                    except ValueError: # Sometimes it is written exponentially for very small separations, so make it zero.
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

    print "Creating -x periodic image..."
    # Create -x image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        xCoords = atom[3]-(abs(volume[0][0])+abs(volume[0][1]))
        if xCoords >= volume[0][0]-border:
            totalAtoms += 1
            imageAtom = [totalAtoms, atom[1], totalChains, xCoords, atom[4], atom[5], atom[6], atom[7], atom[8]]
            xMinusImage.append(imageAtom)
    print "Done!"
    print "Creating +x periodic image..."

    # Create +x image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        xCoords = atom[3]+(abs(volume[0][0])+abs(volume[0][1]))
        if xCoords <= volume[0][1]+border:
            totalAtoms += 1
            imageAtom = [totalAtoms, atom[1], totalChains, xCoords, atom[4], atom[5], atom[6], atom[7], atom[8]]
            xPlusImage.append(imageAtom)
    print "Done!"
    print "Concatenating x images..."
    trajectory = xMinusImage+trajectory+xPlusImage
    print "X-axis periodic images complete!"
    print "Creating -y periodic image..."

    # Create -y image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        yCoords = atom[4]-(abs(volume[1][0])+abs(volume[1][1]))
        if yCoords >= volume[1][0]-border:
            totalAtoms += 1
            imageAtom = [totalAtoms, atom[1], totalChains, atom[3], yCoords, atom[5], atom[6], atom[7], atom[8]]
            yMinusImage.append(imageAtom)
    print "Done!"
    print "Creating +y periodic image..."

    # Create +y image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        yCoords = atom[4]+(abs(volume[1][0])+abs(volume[1][1]))
        if yCoords <= volume[1][1]+border:
            totalAtoms += 1
            imageAtom = [totalAtoms, atom[1], totalChains, atom[3], yCoords, atom[5], atom[6], atom[7], atom[8]]
            yPlusImage.append(imageAtom)
    print "Done!"
    print "Concatenating y images..."
    trajectory = yMinusImage+trajectory+yPlusImage
    print "X and Y-axis periodic images complete!"
    print "Creating -z periodic image..."

    # Create -z image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        zCoords = atom[5]-(abs(volume[2][0])+abs(volume[2][1]))
        if zCoords >= volume[2][0]-border:
            totalAtoms += 1
            imageAtom = [totalAtoms, atom[1], totalChains, atom[3], atom[4], zCoords, atom[6], atom[7], atom[8]]
            zMinusImage.append(imageAtom)
    print "Done!"
    print "Creating +z periodic image..."

    # Create +z image
    for atom in trajectory:
        if atom[2] != previousChainNumber:
            previousChainNumber = atom[2]
            totalChains += 1
        zCoords = atom[5]+(abs(volume[2][0])+abs(volume[2][1]))
        if zCoords <= volume[2][1]+border:
            totalAtoms += 1
            imageAtom = [totalAtoms, atom[1], totalChains, atom[3], atom[4], zCoords, atom[6], atom[7], atom[8]]
            zPlusImage.append(imageAtom)
    print "Done!"
    print "Concatenating z images..."

    trajectory = zMinusImage+trajectory+zPlusImage
    print "All 8 periodic images complete!"
    print "There are now", totalChains, "chains with", totalAtoms, "atoms in the second-order system."

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



def calculateParacrystallinity(atoms, bonds, angles, dihedrals, impropers, trajectory, simVolData, hkl, periodicityOrder, tolerance, crystalCutOff, border, limitThiophenesExamined, thiopheneNumber, blockAdjacent):
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


    # plotGraph(trajectory, 'original')
    # plotGraph(secondOrderPeriodic, 'border')


    for atom in secondOrderPeriodic:
        if atom[1] == 1:
            secondOrderThioAtoms.append(atom)


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


        # For each conjugated subunit, find the pi stacking direction
        # print "EXAMINING THIO ATOM", improper[0]
        # print "Improper =", improper
        
#         angleAtoms = []
#         for angle in thioThioThioAngle:
# #            print angle
#             if angle[1] == thioAtom[0]:
#                 backBoneAngle = angle
#                 # Should only be one of these angles
#                 break
#         try:
#             backBoneAngle = backBoneAngle[:]
# #            print "ANGLE IT BELONGS TO =", backBoneAngle
#         except:
#             # This is a starting or ending monomer so has no line of three thios, also exclude these from the averaging process
#             continue
#         print "P1-P1-P1 Backbone Angle =", backBoneAngle
#         # Get the angle atoms
#         for atomIndex in backBoneAngle:
#             angleAtoms.append(thioAtoms[((atomIndex+2)/3)-1])
#         print angleAtoms


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


        closestAtom = findClosestAtom(improper[0], secondOrderThioAtoms, stackingVector, tolerance, bondedAtoms)
        if (closestAtom == 0) or (closestAtom[1] >= crystalCutOff):
            thiophenesWithoutNeighbours += 1
            closestAtom = [0, 100, 0, [0,0,0], [0,0,0]]

#        print "Examining Atom", improper[0], "..."


    #     closestAtom = findClosestAtom(improper[0], thioAtoms, normalVector, tolerance)

    #     if (closestAtom == 0):
    #         if (periodicityOrder > 1):
    #             # No nearby atom in this volume - increase the volume to include surrounding periods
    #             if (secondOrderPeriodic == False) and (secondOrderPeriodicCubes == False):
    # #                print "No neighbours in first-order volume, generating second-order periods..."
    #                 secondOrderPeriodic, secondOrderPeriodicCubes, secondOrderVol = generatePeriodicity(trajectory, simVolData)
    #                 for atom in secondOrderPeriodicCubes:
    #                     if atom[1] == 1:
    #                         secondOrderThioAtoms.append(atom)
    #             secondOrderClosestAtom = findClosestAtom(improper[0], secondOrderThioAtoms, normalVector, tolerance)
    #             if (secondOrderClosestAtom == 0):
    #                 if (periodicityOrder > 2):
    #                     if (thirdOrderPeriodic == False) and (thirdOrderPeriodicCubes == False):
    #     #                    print "No neighbours in second-order volume, generating third-order periods..."
    #                         thirdOrderPeriodic, thirdOrderPeriodicCubes, thirdOrderVol = generatePeriodicity(secondOrderPeriodic, secondOrderVol)
    #                         for atom in thirdOrderPeriodicCubes:
    #                             if atom[1] == 1:
    #                                 thirdOrderThioAtoms.append(atom)
    #                     thirdOrderClosestAtom = findClosestAtom(improper[0], thirdOrderThioAtoms, normalVector, tolerance)
    #                     if thirdOrderClosestAtom == 0:
    #                         noNeighbourSeparation = np.sqrt(2*(thirdOrderVol[0][1]**2))
    #                         # print "No third-order neighbours found for atom", str(thioAtom)+"."
    #                         # print "Artificially setting separation to the extent of the third-order periodic boundaries =", noNeighbourSeparation
    #                         closestAtom = [0, noNeighbourSeparation, 0, [0, 0, 0]]
    #                         thiophenesWithoutNeighbours += 1
    #                     else:
    #     #                    print "Third order atom found"
    #                         closestAtom = thirdOrderClosestAtom
    #                 else:
    #                     noNeighbourSeparation = np.sqrt(2*(secondOrderVol[0][1]**2))
    #                     closestAtom = [0, noNeighbourSeparation, 0, [0, 0, 0]]
    #                     thiophenesWithoutNeighbours += 1
    #             else:
    #                 # print "Second order atom found"
    #                 closestAtom = secondOrderClosestAtom
    #         else:
    #             noNeighbourSeparation = np.sqrt(2*(simVolData[0][1]**2))
    #             closestAtom = [0, noNeighbourSeparation, 0, [0, 0, 0]]
    #             thiophenesWithoutNeighbours += 1

#        print closestAtom

        dhkl.append(closestAtom[1])
        dhkl2.append(closestAtom[1]**2)


        atomsTreated += 1


#         if closestAtom[0] != 0:
#             plotThisOne += 1
# #            plotGraph(trajectory, 'test', atom1ToHighlight=improper[0], atom2ToHighlight=closestAtom)
#             if plotThisOne == 10:
#                 plotGraph2(trajectory, improper, closestAtom, normalVector)
#                 raise SystemError('Custard')




        # if len(dhkl) != atomsTreated:
        #     print "Len(dhkl) =", len(dhkl)
        #     print "Atoms treated =", atomsTreated
        #     print "Thio Atom =", thioAtom
        #     print "Backbone Angle =", backBoneAngle
        #     print "Normal Vector =", normalVector
        #     raw_input('DHKL length not equal to the number of atoms treated!! Press return to continue...')

    # Now we have dhkl for all of the atoms, let's take an average for this timestep
#    avdhkl2 = np.sum(dhkl2)/len(dhkl2)#atomsTreated
#    avdhkl = np.sum(dhkl)/len(dhkl)#atomsTreated

#    avdhkl = np.mean(dhkl)
#    avdhkl2 = np.mean(dhkl2)

    avdhkl = np.sum(dhkl)/atomsTreated
    avdhkl2 = np.sum(dhkl2)/atomsTreated

    g2 = (avdhkl2 - avdhkl**2)/(avdhkl**2)


    print "\n-----=====================-----"
    print "Paracrystallinity parameter g =", np.sqrt(g2)
    print "Proportion of thiophenes with no neighbours (calcs truncated) =", float(thiophenesWithoutNeighbours)/float(atomsTreated)
    print "-----=====================-----\n"

    # P.figure()
    # x = np.arange(0, int(max(dhkl)))
    # P.hist(dhkl, x)
    # P.show()

    # P.figure()
    # x = np.arange(0, int(max(dhkl2)))
    # P.hist(dhkl2, x)
    # P.show()

    return np.sqrt(g2), dhkl

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



def findClosestAtom(thioAtom, allThioAtoms, axis, toleranceAngle, bondedAtoms): # ANGULAR
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


### THIS IS ALL WRONG SO IGNORE IT COMPLETELY
# #         # Do this by angles instead. Cross product of separation vector and axis vector divided by the magnitudes of both
# #         # == 1 = sin(theta).
# #         # The closer this value is to zero, the closer the two vectors are, so sort by theta and take the first one.
# #         # If theta <= 90 then theta = theta
# #         # If 90 < theta <= 180 then theta = 180-theta (close to plane on reverse side of the molecule)
# #         # If 180 < theta <= 270 then theta = theta-180 (behind molecule, reverse side)
# #         # If theta > 270 then theta = 360-theta (behind molecule)
# #         # D then becomes its separation.

#         crossProductMagnitude = np.linalg.norm(np.cross(separationVector, axis))
# #        print "\nnorm =", crossProductMagnitude
# #        print crossProductMagnitude - 1
#         # For some strange reason, the program outputs infinity if you do arcsin(1), even though the individual steps
#         # work fine and output the correct numbers....put in a catch for that:
#         if abs(crossProductMagnitude - 1) < 1E-15:
#             vectorAngle = np.pi/2.
#         else:
#             vectorAngle = np.arcsin(crossProductMagnitude)


#         vectorAngleRaw = vectorAngle

#         if vectorAngle <= (np.pi/2.):
#             vectorAngle = vectorAngle
#         elif vectorAngle <= (np.pi):
#             vectorAngle = (np.pi/2.) - vectorAngle
#         elif vectorAngle <= (3*np.pi/2.):
#             vectorAngle = vectorAngle - (np.pi/2.)
#         elif vectorAngle <= (2*np.pi):
#             vectorAngle = np.pi - vectorAngle
#         else:
#             print "---======---"
#             print "Atom =", thioAtom
#             print "Atom2 =", atom2
#             print "Separation Vector =", separationVector
#             print "Axis =", axis
#             print "Cross product =", np.cross(separationVector, axis)
#             print "Normalised Cross Product =", np.linalg.norm(np.cross(separationVector, axis))
#             print "VectorAngle =", vectorAngle
#             raise SystemError('Angle calculation is incorrect')


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
            closestAtomsList.append([atom2Index, calculateSeparation(atom1Position, atom2Position), dotProductAngle, separationVector, atom2Position])

    closestAtomsList.sort(key = lambda x: x[1])
#    print closestAtomsList

#    print "Length =", closestAtomsList

    # ClosestAtomsList is in the format: [ Atom Index, Separation, AngleBetweenSeparationVectorAndAxis, [SepVecX, SepVecY, SepVecZ] ]
    if len(closestAtomsList) != 0:
        return closestAtomsList[0]
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
    

def generatePeriods(trajectory, originalVolume, periodsToInclude):
    totalAtoms = trajectory[-1][0]
    totalChains = trajectory[-1][2]
    previousChainNumber = trajectory[-1][2]
    newVolume = []
    xMinusImage = []
    xPlusImage = []
    yMinusImage = []
    yPlusImage = []
    zMinusImage = []
    zPlusImage = []
    newVolume = originalVolume[:]
    # Need to include all possible combinations of permitted periods!
    while np.sum(periodsToInclude) >= 1:
        if periodsToInclude[0] == 1:
            # Do x-minus
            print "Creating x-minus period..."
            for atom in trajectory:
                totalAtoms += 1
                if atom[2] != previousChainNumber:
                    previousChainNumber = atom[2]
                    totalChains += 1
                imageAtom = [totalAtoms, atom[1], totalChains, atom[3]-(abs(newVolume[0][0])+abs(newVolume[0][1])), atom[4], atom[5], atom[6], atom[7], atom[8]]
                xMinusImage.append(imageAtom)
            trajectory += xMinusImage
            periodsToInclude[0] = 0
            newVolume[0][0] -= (abs(newVolume[0][0])+abs(newVolume[0][1]))
        if periodsToInclude[1] == 1:
            # Do x-plus
            print "Creating x-plus period..."
            for atom in trajectory:
                totalAtoms += 1
                if atom[2] != previousChainNumber:
                    previousChainNumber = atom[2]
                    totalChains += 1
                imageAtom = [totalAtoms, atom[1], totalChains, atom[3]+(abs(newVolume[0][0])+abs(newVolume[0][1])), atom[4], atom[5], atom[6], atom[7], atom[8]]
                xPlusImage.append(imageAtom)
            trajectory += xPlusImage
            periodsToInclude[1] = 0
            newVolume[0][1] += (abs(newVolume[0][0])+abs(newVolume[0][1]))
        if periodsToInclude[2] == 1:
            # Do y-minus
            print "Creating y-minus period..."
            for atom in trajectory:
                totalAtoms += 1
                if atom[2] != previousChainNumber:
                    previousChainNumber = atom[2]
                    totalChains += 1
                imageAtom = [totalAtoms, atom[1], totalChains, atom[3], atom[4]-(abs(newVolume[0][0])+abs(newVolume[0][1])), atom[5], atom[6], atom[7], atom[8]]
                yMinusImage.append(imageAtom)
            trajectory += yMinusImage
            periodsToInclude[2] = 0
            newVolume[1][0] -= (abs(newVolume[1][0])+abs(newVolume[1][1]))
        if periodsToInclude[3] == 1:
            # Do y-plus
            print "Creating y-plus period..."
            for atom in trajectory:
                totalAtoms += 1
                if atom[2] != previousChainNumber:
                    previousChainNumber = atom[2]
                    totalChains += 1
                imageAtom = [totalAtoms, atom[1], totalChains, atom[3], atom[4]+(abs(newVolume[0][0])+abs(newVolume[0][1])), atom[5], atom[6], atom[7], atom[8]]
                yPlusImage.append(imageAtom)
            trajectory += yPlusImage
            periodsToInclude[3] = 0
            newVolume[1][1] += (abs(newVolume[1][0])+abs(newVolume[1][1]))
        if periodsToInclude[4] == 1:
            # Do z-minus
            print "Creating z-minus period..."
            for atom in trajectory:
                totalAtoms += 1
                if atom[2] != previousChainNumber:
                    previousChainNumber = atom[2]
                    totalChains += 1
                imageAtom = [totalAtoms, atom[1], totalChains, atom[3], atom[4], atom[5]-(abs(newVolume[0][0])+abs(newVolume[0][1])), atom[6], atom[7], atom[8]]
                zMinusImage.append(imageAtom)
            trajectory += zMinusImage
            periodsToInclude[4] = 0
            newVolume[2][0] -= (abs(newVolume[2][0])+abs(newVolume[2][1]))
        if periodsToInclude[5] == 1:
            # Do z-plus
            print "Creating z-plus period..."
            for atom in trajectory:
                totalAtoms += 1
                if atom[2] != previousChainNumber:
                    previousChainNumber = atom[2]
                    totalChains += 1
                imageAtom = [totalAtoms, atom[1], totalChains, atom[3], atom[4], atom[5]+(abs(newVolume[0][0])+abs(newVolume[0][1])), atom[6], atom[7], atom[8]]
                zPlusImage.append(imageAtom)
            trajectory += zPlusImage
            periodsToInclude[5] = 0
            newVolume[2][1] += (abs(newVolume[2][0])+abs(newVolume[2][1]))

    print "Original Volume =", originalVolume
    print "New Volume =", newVolume

    return trajectory, newVolume
        


def findIndex(string, logical):
    '''This function returns the locations of an inputted character (logical) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == logical:
            locations.append(index)
        index += 1
    return locations
    
def writeCSV(data, name):
    '''Writes the dhkl values in a CSV file'''
    filename = './outputGraphs/'+name+'.csv'
    document = csv.writer(open(filename, 'w+'), delimiter = ',')
    for dhkl in data:
        document.writerow([dhkl])
    print 'Dhkl data written to', filename


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
    parser.add_argument("-hkl", "--hkl", default="[0, 1, 0]", help="The hkl direction to calculate the paracrystallinity for. [0, 1, 0] is pi-stacking direction, [0, 0, 1] is along-chain direction.")
    args = parser.parse_args()


#    dumpstepsToBeTreated = np.arange(0, 501, 50)
    # Every dumpstep is at 2,000 lammps timesteps
#    dumpstepsToBeTreated = [0]

    # dumpstepsToBeTreated = np.arange(0, 501, 50)
    paraData = []
    paraDataSmallDist = []
    periodicityOrder = 3 # The maximum number of periodic orders to consider before setting dhkl to an artificial maximum based on the simulation volume

    trajectoryNameRaw = args.lammpstrj
    datNameRaw = args.dat
    dumpstep = int(args.dumpstep)
    hkl = []
    for axisValue in args.hkl.split():
        valueOnly = ""
        for character in axisValue:
            try:
                int(character)
                valueOnly += character
            except:
                continue
        hkl.append(int(valueOnly))

    if (hkl == [0, 1, 0]):
        hklString = 'pi'
    elif (hkl == [0, 0, 1]):
        hklString = 'ch'
    elif (hkl == [1, 0, 0]):
        hklString = 'al'
    else:
        hklString = 'xx'

    direc = './'
    tolerance = np.pi/6.
    crystalCutOff = 999999
    border = 20
    limitThiophenesExamined = False
    blockAdjacent = True
    thiopheneNumber = 100

    print "Tolerance =", tolerance

    trajectoryName = direc+trajectoryNameRaw
    datName = direc+datNameRaw


    trajectoryData, allTimestepNos, simVolData = loadTrajectory(trajectoryName)
    masses, atoms, bonds, angles, dihedrals, impropers = loadDat(datName)
    timestepNos = []

    if dumpstep >= len(trajectoryData):
        # Just run the final dumpstep
        dumpstep = -1


    print "Treating dumpstep", str(dumpstep)+"..."
    # periodicTraj = generatePeriodicity(trajectoryData[dumpstep]. simVolData[dumpstep])
    paraG, dhkl = calculateParacrystallinity(atoms, bonds, angles, dihedrals, impropers, trajectoryData[dumpstep], simVolData[dumpstep], hkl, periodicityOrder, tolerance, crystalCutOff, border, limitThiophenesExamined, thiopheneNumber, blockAdjacent)
    paraData.append(paraG)
    timestepNos.append(allTimestepNos[dumpstep])
    dumpNumber = str(dumpstep+1)
    if dumpNumber == '0':
        dumpNumber = str(len(trajectoryData))
    while len(dumpNumber) < 3:
        dumpNumber = '0'+dumpNumber
    
    name = datNameRaw[:-4]+'_step'+dumpNumber+'_hkl_'+hklString
    writeCSV(dhkl, name)
    paraGSmallDist = plotHist(dhkl, name)
    paraDataSmallDist.append(paraGSmallDist)

    dotsList = findIndex(trajectoryNameRaw, '.')
    graphName = trajectoryNameRaw[:dotsList[1]]

    plot(timestepNos, paraData, graphName)
    plot(timestepNos, paraDataSmallDist, graphName+'_smallDist')







#     if len(dumpstepsToBeTreated) == 0:
#         for trajNo in range(len(trajectoryData)):
#             print "Treating dumpstep", trajNo
# #            periodicTraj = generatePeriodicity(trajectoryData[trajNo], simVolData[trajNo])
#             paraG, dhkl = calculateParacrystallinity(atoms, bonds, angles, dihedrals, impropers, trajectoryData[trajNo], simVolData[trajNo], periodicityOrder, tolerance, crystalCutOff, border, limitThiophenesExamined, thiopheneNumber, blockAdjacent)
#             paraData.append(paraG)
#             timestepNos.append(allTimestepNos[trajNo])
#     else:
#         for dumpstep in dumpstepsToBeTreated:
