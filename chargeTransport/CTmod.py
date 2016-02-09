import os
import numpy as np
import random as R
import pylab as P
import mpl_toolkits.mplot3d.axes3d as p3
import copy
import time as T
import csv
import subprocess
import argparse
import sys



elementaryCharge = 1.60217657E-19 # C
kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34 # m^{2} kg s^{-1}


def loadDat(morphology):
    originalFile = open('./'+str(morphology)+'.dat', 'r')
    datData = originalFile.readlines()
    originalFile.close()

    simVolRawData = datData[14:17]
    simVolData = []

    for element in simVolRawData:
        temp = element.split(' ')
        simVolData.append([float(temp[0]), float(temp[1])])

    return simVolData

def loadSegmentData(morphology, koopmansApproximation):
    # We want there to be a dictionary, sorted by Host Segments that looks like the following:
    # {HostSegmentNo: [Hop Option Numbers, Hop Option Transfer Integral, [Hop option Coords], [Periodic Image Coords]]}}

    unmappedHOMODict = getHOMOData(morphology)                                                        # The raw data from the ZINDO calculations (sigma too large)
    HOMODict = getMappedHOMOData(unmappedHOMODict)                                                    # Of the form {(SegmentID): HOMO}

    # print "TEST:"
    testMean, testSigma = getGaussianStatistics(HOMODict.values())
    # print "testMean =", testMean
    # print "testSigma =", testSigma
    transferIntegralDictionary = getTransferIntegralData(morphology, koopmansApproximation, HOMODict) # Of the form {(Segment1, Segment2): TransferIntegral}
    COMPairDictionary = getCOMPairData(morphology)                                                    # Of the form {(Segment1, Segment2): [COM1, COM2]}
    COMSingleDictionary = getCOMSingleData(morphology)                                                # Of the form {(SegmentID): [PhysicalCOM]}
    segmentLength = getSegmentLength(morphology)                                                      # Of the form {(SegmentID): length}
    totalNumberOfSegments = len(COMSingleDictionary.items())

    hoppingMaster = {}

    for segmentPair in transferIntegralDictionary.keys():
        # Add forward hop first
        forwardHopData = []
        [COM1, COM2] = COMPairDictionary[segmentPair]
        transferIntegral = transferIntegralDictionary[segmentPair]
        # Check if this is a periodic hop!
        # By crossreferencing the COM2 of the hop target segment with the COMSingleDictionary, we can work out if we've just hopped across a boundary
        periodicImageCoords = [0, 0, 0]
        realSegmentCoords = COMSingleDictionary[segmentPair[1]]
        if COM2 != realSegmentCoords:
            # Periodic hop, destination is an image segment
            for coordinate in range(len(COM2)):
                if COM2[coordinate] > realSegmentCoords[coordinate]:
                    # /*It has been moved in the POSITIVE X,Y or Z direction making it in the NEGATIVE X,Y or Z image*/
                    # NB: APPARENTLY THIS IS WRONG AS IT ENDS UP REVERSING THE CORRECT IMAGE COORDINATES
                    periodicImageCoords[coordinate] += 1
                elif COM2[coordinate] < realSegmentCoords[coordinate]:
                    # /*It has been moved in the NEGATIVE X,Y or Z direction making it in the POSITIVE X,Y or Z image*/
                    # NB: APPARENTLY THIS IS WRONG AS IT ENDS UP REVERSING THE CORRECT IMAGE COORDINATES
                    periodicImageCoords[coordinate] -= 1
        forwardHopData = [segmentPair[1], transferIntegral, COM2, periodicImageCoords]
        # Add this to the hoppingMaster
        if (segmentPair[0] in hoppingMaster):
            hoppingMaster[segmentPair[0]].append(forwardHopData)
        else:
            hoppingMaster[segmentPair[0]] = [forwardHopData]
        # Now add the reverse hop
        reverseHopData = [segmentPair[0], transferIntegral, COM1, [-1*periodicImageCoords[0], -1*periodicImageCoords[1], -1*periodicImageCoords[2]]] # Have to reverse the periodic image for when we keep track of the periodic boundary passes
        if (segmentPair[1] in hoppingMaster):
            hoppingMaster[segmentPair[1]].append(reverseHopData)
        else:
            hoppingMaster[segmentPair[1]] = [reverseHopData]

    return hoppingMaster, COMSingleDictionary, segmentLength, HOMODict


def getHOMOData(morphology):
    singleSegmentsCSV = open('./preparedMorphologies/'+morphology+'/SingleSegments.csv', 'r')
    singleSegmentsData = csv.reader(singleSegmentsCSV, delimiter = ',')
    HOMODict = {}
    for row in singleSegmentsData:
        segmentNumber = int(row[0])
        HOMO = float(row[1])
        HOMODict[segmentNumber] = HOMO
    singleSegmentsCSV.close()
    return HOMODict

def getMappedHOMOData(HOMODict):
    mappedHOMODict = {}
    targetSigma = 0.07 # 70 meV
    HOMOMean, HOMOSigma = getGaussianStatistics(HOMODict.values())
    for segmentNumber in HOMODict.keys():
        unmappedHOMO = HOMODict[segmentNumber]
        unmappedHOMOSigmasFromMean = float(abs(HOMOMean - unmappedHOMO))/float(HOMOSigma)
        if unmappedHOMO <= HOMOMean:
            mappedHOMO = HOMOMean - (unmappedHOMOSigmasFromMean * targetSigma)
        else:
            mappedHOMO = HOMOMean + (unmappedHOMOSigmasFromMean * targetSigma)
        mappedHOMODict[segmentNumber] = mappedHOMO
    return mappedHOMODict


def getGaussianStatistics(distribution):
    mean = np.mean(distribution)
    sigma = np.std(distribution)
    return mean, sigma


def getTransferIntegralData(morphology, koopmansApproximation, HOMODict):
    transferIntegralCSV = open('./preparedMorphologies/'+morphology+'/TransferIntegrals.csv', 'r')
    transferIntegralData = csv.reader(transferIntegralCSV, delimiter = ',')
    transferIntegralDict = {}

    numberOfTransferIntegrals = 0
    numberOfKoopmansFails = 0

    for row in transferIntegralData:
        segmentPair = (int(row[0]), int(row[1]))
        numberOfTransferIntegrals += 1
        if koopmansApproximation == True:
            transferIntegral = float(row[4])
        else:
            koopmanTransferIntegral = float(row[4]) # Koopman approximated transfer integral
            HOMOSplitting = 2*koopmanTransferIntegral
            # GET the site Energy difference squared
            energyDifference = (HOMODict[int(row[0])] - HOMODict[int(row[1])])
            if ((HOMOSplitting**2) < (energyDifference**2)):
                numberOfKoopmansFails += 1
                transferIntegralDict[segmentPair] = 0
                continue
            transferIntegral = 0.5*np.sqrt((HOMOSplitting**2) - (energyDifference**2))
        transferIntegralDict[segmentPair] = transferIntegral
    transferIntegralCSV.close()
    if numberOfKoopmansFails > 0:
        print "There are", numberOfTransferIntegrals, "transfer integrals in total, and there are", numberOfKoopmansFails, "koopman-related failures."
    return transferIntegralDict


def getCOMPairData(morphology):
    COMPairCSV = open('./preparedMorphologies/'+morphology+'/SegmentPairCOM.csv', 'r')
    COMPairData = csv.reader(COMPairCSV, delimiter = ',')
    COMPairDict = {}
    for row in COMPairData:
        segmentPair = (int(row[0]), int(row[1]))
        COM1 = [float(row[2]), float(row[3]), float(row[4])]
        COM2 = [float(row[5]), float(row[6]), float(row[7])]
        COMPairDict[segmentPair] = [COM1, COM2]
    COMPairCSV.close()
    return COMPairDict


def getCOMSingleData(morphology):
    COMSingleCSV = open('./preparedMorphologies/'+morphology+'/SegmentSingleCOM.csv', 'r')
    COMSingleData = csv.reader(COMSingleCSV, delimiter = ',')
    COMSingleDict = {}
    for row in COMSingleData:
        segmentID = int(row[0])
        COM = [float(row[1]), float(row[2]), float(row[3])]
        COMSingleDict[segmentID] = COM
    COMSingleCSV.close()
    return COMSingleDict
    


def getSegmentLength(morphology):
    segmentLengthCSV = open('./preparedMorphologies/'+morphology+'/SingleSegments.csv', 'r')
    segmentLengthData = csv.reader(segmentLengthCSV, delimiter = ',')
    segmentLengthDict = {}
    for row in segmentLengthData:
        segmentID = int(row[0])
        length = int(row[2])
        segmentLengthDict[segmentID] = length
    segmentLengthCSV.close()
    return segmentLengthDict



def sortTrajByMolecule(trajectoryData):
    # trajectoryData has the form [atom1, atom2, atom3 .... atomN]
    # where atom has the form [index, type, molecule, x, y, z, ix, iy, iz]

    # sort atoms by molecule
    trajectoryData.sort(key = lambda x: x[2])

    thioAtomsByMolecule = [{}]
    alk1AtomsByMolecule = [{}]
    alk2AtomsByMolecule = [{}]


    moleculeNumber = 1
    
    for atom in trajectoryData:
        if atom[2] != moleculeNumber: # New molecule so add infrastructure to the molecule lists
            thioAtomsByMolecule.append({})
            alk1AtomsByMolecule.append({})
            alk2AtomsByMolecule.append({})
            moleculeNumber = atom[2]
        if atom[1] == 1: # Thio Atom
             thioAtomsByMolecule[atom[2]-1][atom[0]] = [atom[3], atom[4], atom[5]]
#            thioAtomsByMolecule[atom[2]-1].append([atom[0], [atom[3], atom[4], atom[5]]])
        elif atom[1] == 2: # Alk1 Atom
            alk1AtomsByMolecule[atom[2]-1][atom[0]] = [atom[3], atom[4], atom[5]]
#            alk1AtomsByMolecule[atom[2]-1].append([atom[0], [atom[3], atom[4], atom[5]]])
        elif atom[1] == 3: # Alk2 Atom
            alk2AtomsByMolecule[atom[2]-1][atom[0]] = [atom[3], atom[4], atom[5]]
#            alk2AtomsByMolecule[atom[2]-1].append([atom[0], [atom[3], atom[4], atom[5]]])
        else:
            raise SystemError('Unknown atom type "'+str(atom[1])+'" for atom '+str(atom))
        # Atom coordinates are now split into 3 lists depending on their coarse-grain site type in the example form:
        # [thiosInMolecule1, thiosInMolecule2, thiosInMolecule3...]
        # where thiosInMolecule1 has the form {thio1Index: [thio1X, thio1Y, thio1Z], thio2Index: [thio2X, thio2Y, thio2Z],...}


    return thioAtomsByMolecule, alk1AtomsByMolecule, alk2AtomsByMolecule

def randomPosition(simVolData, COMs):
    randX = R.uniform(simVolData[0][0], simVolData[0][1])
    randY = R.uniform(simVolData[1][0], simVolData[1][1]) 
    randZ = R.uniform(simVolData[2][0], simVolData[2][1])
    randomPosn = [randX, randY, randZ]

    separationToSegments = [] # Of the form [ [seg1No, seg1Sep], [seg2No, seg2Sep] ... ]
    for segNo, COM in COMs.items():
        separationToSegments.append([segNo, calcSeparation(randomPosn, COM)])

    separationToSegments.sort(key = lambda x: x[1])
    return int(separationToSegments[0][0])

def determineHopTime(rate):
    # Determine the wait times for a hop using the MC equation based on a given rate
    if rate != 0:
        while True:
            x = R.random()
            print "Random number =", x
            if (x != 0.0) and (x != 1.0):
                break
        tau = -np.log(x)/rate
    else:
        # Rate is zero, therefore infinite time until hop occurs - set it to very large
        tau = 1E20
    return tau


def calcSeparation(coords1, coords2):
    # Calculates the physical separation between two positions of the form
    # [x, y, z]
    separation = np.sqrt((coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 + (coords1[2]-coords2[2])**2)
    return separation



def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    return locations


def writeCSVFile(fileName, carrierNo, displacement):
    fileHandle = open(fileName, 'a+')
    document = csv.writer(fileHandle, delimiter = ',')
    displacementInM = float(displacement)*1E-10 # Convert from angstroms to metres
    document.writerow([carrierNo, displacementInM])
    fileHandle.close()


class chargeCarrier:
    def __init__(self, initialSegment, hoppingMaster, COMSingleDictionary, HOMODictionary, segmentLengths, simVolData, temperature, plottingSubroutines, koopmansApproximation, simulationBoundaries):
        self.initialPosition = COMSingleDictionary[initialSegment]
        self.position = initialSegment
        self.imagePosition = [0, 0, 0]
        self.hoppingMaster = hoppingMaster
        # Remove the COM of this segment from the hopTargets list
        self.COMSingleDictionary = COMSingleDictionary # The locations of the physical segments
        self.HOMO = HOMODictionary
        self.segmentLengths = segmentLengths
        self.temperature = temperature
        self.simVolData = simVolData
        self.plottingSubroutines = plottingSubroutines
        self.koopmansApproximation = koopmansApproximation
        self.globalTime = 0
        self.simulationBoundaries = simulationBoundaries
        self.reinitialise()


    def reinitialise(self):
        self.hopTargets = self.hoppingMaster[self.position]
        #self.hopTargets of the form [segmentNumber, transferIntegral, Coords, PeriodicImage]
        self.COM = self.COMSingleDictionary[self.position]
        self.segmentLength = self.segmentLengths[self.position]
        if self.plottingSubroutines == True:
            plotCarrier(self.COMSingleDictionary, self.position, self.hopTargets, self.imagePosition, self.simVolData, self.initialPosition, self.globalTime, self.simulationBoundaries)


    def calcHopOptions(self):
        # self.hopTargets is now a list of all the possible hopping targets, unsorted
        # For each hop target, need to calculate the rate and therefore the timescale with which a hop can take place
        hopTimes = []
        hopOptionsPosn = {} # Of the form {segmentNo: [COM, image]}
        lambdaij = self.calculateLambdaij()
        for hopTarget in self.hopTargets:
            hopOptionsPosn[hopTarget[0]] = [hopTarget[2], hopTarget[3]]
            transferIntegral = hopTarget[1]*elementaryCharge # in Joules
            if self.koopmansApproximation == False:
                deltaEij = self.calculateEij(hopTarget[0])
            else:
                deltaEij = 0
            hopRate = self.calcHoppingRate(lambdaij, transferIntegral, deltaEij)
            hopTime = determineHopTime(hopRate)
            hopTimes.append([hopTarget[0], hopTime])
            # print "TransferIntegral =", transferIntegral, "Hopping Rate =", hopRate, "Delta Eij =", deltaEij, "Lambdaij =", lambdaij, "Hopping time =", hopTime
        hopTimes.sort(key = lambda x:x[1])
        hopTarget = hopTimes[0]+hopOptionsPosn[hopTimes[0][0]]
        # print "Current Pos =", self.position,
        # print "Target Pos =", hopTarget,
        # print "TI =", transferIntegral,
        # print "Rate =", hopRate,
        # print "Time =", hopTime,
        return hopTarget


    def calculateEij(self, targetSegment):
        Ei = self.HOMO[self.position]
        Ej = self.HOMO[targetSegment]
        EijeV = Ei - Ej
        EijJ = EijeV*elementaryCharge
        return EijJ


    def calculateLambdaij(self):
        # The equation for the internal reorganisation energy was obtained from the data given in
        # Johansson, E and Larsson, S; 2004, Synthetic Metals 144: 183-191.
        # External reorganisation energy obtained from 
        # Liu, T and Cheung, D. L. and Troisi, A; 2011, Phys. Chem. Chem. Phys. 13: 21461-21470
        lambdaExternal = 0.11 # eV
        if self.segmentLength < 12:
            lambdaInternal = 0.20826 - (self.segmentLength*0.01196)
        else:
            lambdaInternal = 0.06474
        lambdaeV = lambdaExternal+lambdaInternal
        lambdaJ = lambdaeV*elementaryCharge
        return lambdaJ



    def calcHoppingRate(self, lambdaij, transferIntegral, deltaEij):
        if self.koopmansApproximation == True:
            kij = (transferIntegral**2/hbar)*np.sqrt((np.pi)/(lambdaij*kB*self.temperature))
        else:
            kij = (transferIntegral**2/hbar)*np.sqrt((np.pi)/(lambdaij*kB*self.temperature))*np.exp(-((deltaEij-lambdaij)**2)/(4*lambdaij*kB*self.temperature))
        return kij


    def determineNextHopType(self):
        if self.nextIntraHopTime < self.nextInterHopTime:
            self.nextHopType = 'Intra'
            self.nextHopTime = self.nextIntraHopTime
        else:
            self.nextHopType = 'Inter'
            self.nextHopTime = self.nextInterHopTime


    def makeHop(self):
        hopTarget = self.calcHopOptions()
        # hopTarget has the form: [segmentNumber, timeToHop, COMPosn, imageLocation]
        hopDistance = np.sqrt(((self.COM[0] - hopTarget[2][0])**2) + ((self.COM[1] - hopTarget[2][1])**2) + ((self.COM[2] - hopTarget[2][2])**2))
        # print "Hop complete. Hopping time =", hopTarget[1], "HoppingDistance =", hopDistance
        # raw_input("Press RETURN to continue...")


        # Update Timings
        self.globalTime += hopTarget[1]
        # print "\n---======---"
        # print "Hopping from", self.position, "at", self.COM, self.imagePosition, "to", hopTarget[0], "at", hopTarget[2], hopTarget[3]
        self.position = hopTarget[0]
        # Update image position
        # For some reason this is the wrong way round. If the particle hops off to the right along the x axis,
        # it appears on the left hand side of the (-1, 0, 0) image. Easy fix is to subtract the hopTarget image to put it in the (1, 0, 0) image instead
        self.imagePosition[0] += hopTarget[3][0]
        self.imagePosition[1] += hopTarget[3][1]
        self.imagePosition[2] += hopTarget[3][2]
        # Update position, segmentLength, plot image
        self.reinitialise()
        # print "Current image Position =", self.imagePosition
        # print "---======---"
        return self.globalTime

        
    def findBondedAtoms(self):
        bondedThiophenes = []
        for bond in self.thioThioBonds:
            if self.currentAtomIndex == bond[0]:
                bondedThiophenes.append(bond[1])
            elif self.currentAtomIndex == bond[1]:
                bondedThiophenes.append(bond[0])
        return bondedThiophenes



def plotCarrier(COMDictionary, currentSegment, hopTargets, imagePosition, simVolData, initialPosition, globalTime, simulationBoundaries):
    xExtent = simulationBoundaries[0]
    yExtent = simulationBoundaries[1]
    zExtent = simulationBoundaries[2]
    neighbouringSegmentCoordsPhysical = []
    neighbouringSegmentCoordsImage = []
    currentSegmentCoords = COMDictionary[currentSegment]
    for hopTarget in hopTargets:
        segmentNo = hopTarget[0]
        if hopTarget[3] == [0, 0, 0]:
            neighbouringSegmentCoordsPhysical.append(COMDictionary[segmentNo])
        else:
            neighbouringSegmentCoordsImage.append(COMDictionary[segmentNo])
    xNeighbourPhys, yNeighbourPhys, zNeighbourPhys = zip(*neighbouringSegmentCoordsPhysical)
    if len(neighbouringSegmentCoordsImage) != 0:
        xNeighbourIm, yNeighbourIm, zNeighbourIm = zip(*neighbouringSegmentCoordsImage)
        xNeighbourIm = list(xNeighbourIm)
        yNeighbourIm = list(yNeighbourIm)
        zNeighbourIm = list(zNeighbourIm)
    xNeighbourPhys = list(xNeighbourPhys)
    yNeighbourPhys = list(yNeighbourPhys)
    zNeighbourPhys = list(zNeighbourPhys)

    fig = P.figure()
    ax = p3.Axes3D(fig)
    for i in range(len(xNeighbourPhys)):
        ax.scatter(xNeighbourPhys[i], yNeighbourPhys[i], zNeighbourPhys[i], s = 20, c = 'g')
    if len(neighbouringSegmentCoordsImage) != 0:
        for i in range(len(xNeighbourIm)):
            ax.scatter(xNeighbourIm[i], yNeighbourIm[i], zNeighbourIm[i], s = 20, c = 'b')

    ax.scatter(currentSegmentCoords[0], currentSegmentCoords[1], currentSegmentCoords[2], s = 50, c = 'r')

    ax.plot([simVolData[0][0], simVolData[0][1]], [simVolData[1][0], simVolData[1][0]], [simVolData[2][0], simVolData[2][0]], c = 'k')
    ax.plot([simVolData[0][0], simVolData[0][1]], [simVolData[1][1], simVolData[1][1]], [simVolData[2][0], simVolData[2][0]], c = 'k')
    ax.plot([simVolData[0][0], simVolData[0][1]], [simVolData[1][0], simVolData[1][0]], [simVolData[2][1], simVolData[2][1]], c = 'k')
    ax.plot([simVolData[0][0], simVolData[0][1]], [simVolData[1][1], simVolData[1][1]], [simVolData[2][1], simVolData[2][1]], c = 'k')

    ax.plot([simVolData[0][0], simVolData[0][0]], [simVolData[1][0], simVolData[1][1]], [simVolData[2][0], simVolData[2][0]], c = 'k')
    ax.plot([simVolData[0][0], simVolData[0][0]], [simVolData[1][0], simVolData[1][1]], [simVolData[2][1], simVolData[2][1]], c = 'k')
    ax.plot([simVolData[0][1], simVolData[0][1]], [simVolData[1][0], simVolData[1][1]], [simVolData[2][0], simVolData[2][0]], c = 'k')
    ax.plot([simVolData[0][1], simVolData[0][1]], [simVolData[1][0], simVolData[1][1]], [simVolData[2][1], simVolData[2][1]], c = 'k')

    ax.plot([simVolData[0][0], simVolData[0][0]], [simVolData[1][0], simVolData[1][0]], [simVolData[2][0], simVolData[2][1]], c = 'k')
    ax.plot([simVolData[0][0], simVolData[0][0]], [simVolData[1][1], simVolData[1][1]], [simVolData[2][0], simVolData[2][1]], c = 'k')
    ax.plot([simVolData[0][1], simVolData[0][1]], [simVolData[1][0], simVolData[1][0]], [simVolData[2][0], simVolData[2][1]], c = 'k')
    ax.plot([simVolData[0][1], simVolData[0][1]], [simVolData[1][1], simVolData[1][1]], [simVolData[2][0], simVolData[2][1]], c = 'k')

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.set_xlim(-60.0, 60.0)
    ax.set_ylim(-60.0, 60.0)
    ax.set_zlim(-60.0, 60.0)


    currentPosition = [currentSegmentCoords[0]+(imagePosition[0]*xExtent), currentSegmentCoords[1]+(imagePosition[1]*yExtent), currentSegmentCoords[2]+(imagePosition[2]*zExtent)]

    displacement = calcSeparation(currentPosition, initialPosition)
    displacement = list(str(displacement))
    displacementString = ''
    for i in range(5):
        try:
            displacementString += displacement[i]
        except IndexError:
            break

    P.title(str(imagePosition)+", Disp = "+str(displacementString)+", t = "+str(globalTime))
    filelist = os.listdir('./')
    filenameCounter = 0
    for files in filelist:
        if ".png" in files:
            filenameCounter += 1
    filenameAddon = str(filenameCounter)
    while len(filenameAddon) < 3:
        filenameAddon = '0'+filenameAddon

    P.savefig('./test'+filenameAddon+'.png')
    print "File saved as ./test"+filenameAddon+".png"



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--morphology", help="The morphology directory in ./preparedMorphologies, containing the following CSV files:\n\nSegmentPairCOM.csv\nSegmentSingleCOM.csv\nSingleSegments.csv\nTransferIntegrals.csv\n\nwhich are all obtained through the runSegments.py -> runOrca.py -> analyseOrca.py system. Also required is a .dat file in the home directory with the same name as the morphology directory, in order to obtain the simulation dimensions (required for the periodic boundary condition treatment).")
    parser.add_argument("-t", "--time", help="The simulation time at which the hopping processes are terminated and the final displacement calculated.")
    parser.add_argument("-c", "--holeNo", help="The identification number of the carrier that is running through the device for the given --time in this --morphology.")
    # parser.add_argument("-n", "--cores", help="The total number of simultaneous carriers to be performed on Hamilton.")
    # parser.add_argument("-q", "--queue", help="The queue to submit the Hamilton jobs to.")
    parser.add_argument("-p", "--plotting", help="Flag determining whether plotting subroutines are used (unrecommended due to wall clock time). 1 == True, 0 == False")
    parser.add_argument("-k", "--koopmans", help="Flag to determine whether Koopmans Approximation is to be used (site energies ignored for the the transfer integral/rate calculation) or not. 1 == True, 0 == False")
    args = parser.parse_args()

    morphology = args.morphology
    simulationTime = float(args.time)
    carrierNo = int(args.holeNo)
    # coresToUse = int(args.cores)
    # queue = args.queue
    koopmansApproximation = bool(int(args.koopmans))
    plotting = bool(int(args.plotting))

    # print "MORPHOLOGY =", morphology
    # print "SIM TIME =", simulationTime


    # Determine the temperature of the morphology
    underscoreLoc = findIndex(morphology, '_')
    temperature = int(morphology[underscoreLoc[-1]+1:-1])

    # Load in the simulation boundary data
    simVolData = loadDat(morphology)
    xExtent = simVolData[0][1] - simVolData[0][0]
    yExtent = simVolData[1][1] - simVolData[1][0]
    zExtent = simVolData[2][1] - simVolData[2][0]

    # Load the information from the CSV files
    hoppingMaster, COMSingleDictionary, segmentLength, HOMODictionary = loadSegmentData(morphology, koopmansApproximation)

    # Pick a random segment to begin the simulation
    initialSegment = randomPosition(simVolData, COMSingleDictionary) # Just an integer of the closest segment to a random posn within the morphology
    initialPosition = COMSingleDictionary[initialSegment]
    currentPosition = initialPosition

    # Initialise a carrier
    hole = chargeCarrier(initialSegment, hoppingMaster, COMSingleDictionary, HOMODictionary, segmentLength, simVolData, temperature, plotting, koopmansApproximation, [xExtent, yExtent, zExtent])
    numberOfHops = 0

    # Hop for the required time
    while True:
        newGlobalTime = hole.makeHop()
        if plotting == True:
            if numberOfHops == 100:
                break
        else:
            if newGlobalTime > simulationTime:
                break
        numberOfHops += 1
        currentPosition = COMSingleDictionary[hole.position]
        imageLocation = hole.imagePosition
        currentPosition = [COMSingleDictionary[hole.position][0]+(hole.imagePosition[0]*xExtent), COMSingleDictionary[hole.position][1]+(hole.imagePosition[1]*yExtent), COMSingleDictionary[hole.position][2]+(hole.imagePosition[2]*zExtent)]

    displacement = calcSeparation(currentPosition, initialPosition)

    # Update CSV file
    if plotting == False:
        csvFileName = "./CTOutput/"+morphology+"_"+str(simulationTime)+"_Koop"+str(koopmansApproximation)[0]+".csv"
        print numberOfHops, "hops complete. Writing displacement of", displacement, "for carrier number", carrierNo, " in", csvFileName
        writeCSVFile(csvFileName, carrierNo, displacement)
    else:
        print numberOfHops, "hops complete. Graphs plotted. Simulation terminating. No CSV data will be saved while plotting == True."
