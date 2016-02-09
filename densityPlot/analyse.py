import os
import numpy as np
import pylab as P
import random as R
import csv
import mpl_toolkits.mplot3d.axes3d as p3
import time as T
import copy

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


def calculateParacrystallinity(atoms, bonds, angles, dihedrals, impropers, trajectories):
    # Thio atom == 1
    # Thio-Thio bond == 1
    # Thio-Thio-Thio angle == 1
    # ATOM TYPE: Molecular
    # [ atomID, moleculeID, atomType, x, y, z ]
    # Other types:
    # [ ID, Type, atom1, atom2, ...]
    thioAtomNos = []
    for atom in atoms:
         if atom[2] == 1:
             thioAtomNos.append(atom[0])

    thioAtoms = []

    thioThioBonds = []
    for bond in bonds:
        if bond[1] == 1:
            thioThioBonds.append([bond[2], bond[3]])
    
    thioThioThioAngle = []
    for angle in angles:
        if angle[1] == 1:
            thioThioThioAngle.append([angle[2], angle[3], angle[4]])

    for atom in trajectories:
        if atom[1] == 1:
            thioAtoms.append(atom)

    atomsToExamine = []
    # atomsToExamine = [atomNumber, Position, Thio-Thio-ThioPlane]


    thioAtoms.sort(key=lambda x: x[0])


#### ALONG CHAIN DIRECTION
    # How best to do this? Surely if we take the direction of the chain at each thiophene (by looking at the neighbouring bond)
    # the next thiophene in the chain will always be the closest? Maybe that's not a problem?
    # We could take the start monomer to the end monomer, but a lot of the chains are extensively twisted so I'm not so sure
    # how useful that information would be.
    # According to Poelking, d_{001} is along the chain direction and varies from 20 to 40 Angstroms which is about right.
        

#### PI STACKING DIRECTION
    dhkl = []
    dhkl2 = []

    time1 = T.time()
    totalAtoms = len(thioAtoms)
    onePercentAtoms = int(totalAtoms/100)
    atomsTreated = 0
    for thioAtom in thioAtoms:
        if np.remainder(atomsTreated, onePercentAtoms) == 0:
            percentageComplete = int(atomsTreated/onePercentAtoms)
            if percentageComplete != 0 and percentageComplete <= 100:
                print "Calculations", percentageComplete, "percent completed..."
            time2 = T.time()
            if percentageComplete != 0 and percentageComplete <= 100:
                totalTime = (((time2-time1)*100)/percentageComplete)
                timeRemaining = totalTime - (time2-time1)
                print "Remaining time estimate: %.1f seconds." % (timeRemaining)


        # For each conjugated subunit, find the pi stacking direction
#        print "EXAMINING THIO ATOM", thioAtom[0]
        angleAtoms = []
        for angle in thioThioThioAngle:
#            print angle
            if angle[1] == thioAtom[0]:
                backBoneAngle = angle
                # Should only be one of these angles
                break
        try:
            backBoneAngle = backBoneAngle[:]
#            print "ANGLE IT BELONGS TO =", backBoneAngle
        except:
            # This is a starting or ending monomer so has no line of three thios, also exclude these from the averaging process
            continue
        # Get the angle atoms
        for atomIndex in backBoneAngle:
            angleAtoms.append(thioAtoms[((atomIndex+2)/3)-1])
        normalVector = piStackingDirection(thioAtom, angleAtoms)
#        print normalVector

        closestAtom = findClosestAtoms(thioAtom, thioAtoms, normalVector)
        if closestAtom != 0:
            dhkl.append(closestAtom[1])
            dhkl2.append(closestAtom[1]**2)
        atomsTreated += 1

    # Now we have dhkl for all of the atoms, let's take an average for this timestep
    avdhkl2 = np.sum(dhkl2)/atomsTreated
    avdhkl = np.sum(dhkl)/atomsTreated

    g2 = (avdhkl2 - avdhkl**2)/(avdhkl**2)

    print "\n-----=====================-----"
    print "Paracrystallinity parameter g =", np.sqrt(g2)
    print "-----=====================-----\n"

    return np.sqrt(g2)

#    print len(atoms)
#    print len(thioAtoms)


def findClosestAtoms(atom, atomList, axis):
    toleranceAngle = np.pi/180 # Monomers subtending an angle of toleranceAngle degrees are still classed as in the same plane.
#    toleranceDistance = 2 # Distance a conjugated subunit can be from the plane to still be classed as in that plane (angstroms)
#    print "AXIS =", axis
    closestAtomsList = []
    atom1Index = atom[0]
    atom1Position = [atom[3], atom[4], atom[5]]
    for atom2 in atomList:
        if atom == atom2:
            continue
        atom2Index = atom2[0]
        atom2Position = [atom2[3], atom2[4], atom2[5]]
        separationVector = np.array([atom2[3] - atom[3], atom2[4] - atom[4], atom2[5] - atom[5]])
        separationVector = list(normaliseVec(separationVector))

#         # Do this by angles instead. Cross product of separation vector and axis vector divided by the magnitudes of both
#         # == 1 = sin(theta).
#         # The closer this value is to zero, the closer the two vectors are, so sort by theta and take the first one.
#         # If theta <= 90 then theta = theta
#         # If 90 < theta <= 180 then theta = 180-theta (close to plane on reverse side of the molecule)
#         # If 180 < theta <= 270 then theta = theta-180 (behind molecule, reverse side)
#         # If theta > 270 then theta = 360-theta (behind molecule)
#         # D then becomes its separation.

        crossProductMagnitude = np.linalg.norm(np.cross(separationVector, axis))
#        print "\nnorm =", crossProductMagnitude
#        print crossProductMagnitude - 1
        # For some strange reason, the program outputs infinity if you do arcsin(1), even though the individual steps
        # work fine and output the correct numbers....put in a catch for that:
        if abs(crossProductMagnitude - 1) < 1E-15:
            vectorAngle = np.pi/2.
        else:
            vectorAngle = np.arcsin(crossProductMagnitude)
        if vectorAngle <= (np.pi/2.):
            vectorAngle = vectorAngle
        elif vectorAngle <= (np.pi):
            vectorAngle = (np.pi/2.) - vectorAngle
        elif vectorAngle <= (3*np.pi/2.):
            vectorAngle = vectorAngle - (np.pi/2.)
        elif vectorAngle <= (2*np.pi):
            vectorAngle = np.pi - vectorAngle
        else:
            print "Atom =", atom
            print "Atom2 =", atom2
            print "Separation Vector =", separationVector
            print "Axis =", axis
            print "Cross product =", np.cross(separationVector, axis)
            print "Normalised Cross Product =", np.linalg.norm(np.cross(separationVector, axis))
            print "VectorAngle =", vectorAngle
            raise SystemError('Angle calculation is incorrect')

        if vectorAngle <= toleranceAngle:
            # This monomer is within the axis plane therefore add it to the closestAtomsList
            closestAtomsList.append([atom2Index, calculateSeparation(atom1Position, atom2Position), vectorAngle, separationVector])


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

def piStackingDirection(atom, angle):
    # Normal direction to a plane defined by points P1, P2 and P3 =
    # (P3 - P1) vectorProduct (P2 - P1)
    # print "Atom =", atom
    # print "Angle =", angle
    vec1 = np.array([(angle[2][3] - angle[0][3]), (angle[2][4] - angle[0][4]), (angle[2][5] - angle[0][5])])
    vec2 = np.array([(angle[1][3] - angle[0][3]), (angle[1][4] - angle[0][4]), (angle[1][5] - angle[0][5])])
    normalVec = np.cross(vec1, vec2)
    normalVec = normaliseVec(normalVec)
    return normalVec
    


def findIndex(string, logical):
    '''This function returns the locations of an inputted character (logical) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == logical:
            locations.append(index)
        index += 1
    return locations
    

def plot(xData, yData, filename, plotFlag):
    P.figure()
    P.plot(xData, yData)
    P.xlabel('Timestep')
    if plotFlag == 'para':
        P.ylabel('Paracrystallinity')
        P.savefig('./para_'+str(filename)+'.png')
        P.show()
    elif plotFlag == 'dens':
        P.ylabel("Density, gcm^{-3}")
        P.savefig('./dens_'+str(filename)+'.png')
        P.show()



if __name__ == '__main__':
    # Include a thing like the runLammps.py where it automatically works out which file to treat
    validDats, validLammps = getFilesList('./')
    filesList = zip(validDats, validLammps)
    
    print "Valid .dat files (that have corresponding lammps trajectories:", filesList

    calcPara = 0
    calcDensity = 1
    plotPara = 0
    plotDensity = 1

    if len(filesList) == 1:
        datName = filesList[0][0]
        trajectoryName = filesList[0][1]
    else:
        while True:
            dataFile = str(raw_input('Please enter the name of the .dat file to be submitted: '))
            for fileNo in range(len(validDats)):
                if dataFile in validDats[fileNo]:
                    dataFile = validDats[fileNo]
                    foundFile = 1
                    break
            if foundFile == 1:
                print "Using datafile:", validDats[fileNo]
                datName = validDats[fileNo]
                trajectoryName = validLammps[fileNo]
                break
            else:
                print ".dat file not found, please try again..."

    # trajectoryName = '0128P3HT_1.dat.lammpstrj'
    # datName = '0128P3HT_1.dat'
    paraData = []
    densityData = []



    trajectoryData, timestepNos, simVolData = loadTrajectory(trajectoryName)
    masses, atoms, bonds, angles, dihedrals, impropers = loadDat(datName)


    # print simVolData[:10]
    # print timestepNos[:10]

    # raise SystemError('CUSTARD')

    if calcPara == 1:
        for traj in trajectoryData:
            paraG = calculateParacrystallinity(atoms, bonds, angles, dihedrals, impropers, traj)
            paraData.append(paraG)

    if calcDensity == 1:
        for trajNo in range(len(trajectoryData)):
            volume = simVolData[trajNo]
            traj = trajectoryData[trajNo]
            density = calculateDensity(masses, traj, volume)
            print "Timestep number", str(trajNo), "has density =", str(density), "gcm^{-3]"
            densityData.append(density)

    if plotPara == 1:
        plot(timestepNos, paraData, datName[:-4], 'para')

    if plotDensity == 1:
        plot(timestepNos, densityData, datName[:-4], 'dens')
