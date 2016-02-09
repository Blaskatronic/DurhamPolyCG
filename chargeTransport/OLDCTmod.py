import os
import numpy as np
import random as R
import pylab as P
import mpl_toolkits.mplot3d.axes3d as p3
import copy


globalTime = 0
plottingSubroutines = True



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



def randomPosition(moleculeSortedAtoms):
    # moleculeSortedAtoms has the form [[{mol1atom1Index, [mol1atom1Coords], mol1atom2Index, [mol1atom2Coords],...}, {mol2atom1Index...}]

    # pick a molecule at random
    molNo = R.randint(0, len(moleculeSortedAtoms)-1) # Randint is inclusive so include 0 but not len(moleculeSortedAtoms)
    molecule = moleculeSortedAtoms[molNo]


    # pick an atom at random within this molecule
    atomNo = R.choice(molecule.keys())

    # Return the coordinates of this atom as the randomPosition
    randomCoords = molecule[atomNo]

    return randomCoords, molecule


def determineHopTime(rate):
    # Determine the wait times for a hop using the MC equation based on a given rate
    x = R.random()
    if rate != 0:
        tau = -np.log(x)/rate
    else:
        # Rate is zero, therefore infinite time until hop occurs - set it to very large
        tau = 1E20
    return tau


class chargeCarrier:
    def __init__(self, initialPosition, thioAtomsByMolecule, molecule, bonds):
        self.position = initialPosition
        self.allThios = thioAtomsByMolecule
        self.thiosInThisMolecule = molecule
        currentAtom = 0
        for atomIndex, atomCoords in molecule.iteritems():
            if atomCoords == initialPosition:
                self.currentAtomIndex = atomIndex
                break


        # Bonds is of the form [bondIndex, bondType, atom1, atom2]
        # where bondtypes are defined in LAMMPS (bondType 1 is Thio-Thio)
        self.thioThioBonds = []
        for bond in bonds:
            if bond[1] == 1: # Thio-thio bond
                self.thioThioBonds.append([bond[2], bond[3]])
        self.intraChainRate = 100.0
        self.interChainRate = 0
        self.nextIntraHopTime = determineHopTime(self.intraChainRate)
        self.nextInterHopTime = determineHopTime(self.interChainRate)
        self.determineNextHopType()
        if plottingSubroutines == True:
            plotCarrier(self.position, self.thiosInThisMolecule, self.currentAtomIndex)


    def determineNextHopType(self):
        if self.nextIntraHopTime < self.nextInterHopTime:
            self.nextHopType = 'Intra'
            self.nextHopTime = self.nextIntraHopTime
        else:
            self.nextHopType = 'Inter'
            self.nextHopTime = self.nextInterHopTime

        
    def makeHop(self):
        print "time until next hop =", self.nextHopTime
        if self.nextHopType == 'Intra':
            # Stay within chain so find bonded atoms and pick a random adjacent thiophene to jump to within this chain
            adjacentThios = self.findBondedAtoms()
            if R.random() <= 0.5:
                hopTarget = adjacentThios[0]
            else:
                hopTarget = adjacentThios[1]
        else:
            pass

        # Hop target is now obtained, so update the various timings and change the position of this charge carrier
        global globalTime
        globalTime += self.nextHopTime
        self.nextIntraHopTime -= self.nextHopTime
        self.nextInterHopTime -= self.nextHopTime

        if self.nextIntraHopTime == 0: # This might have to be changed to <= 1E-20 or something due to rounding errors
            self.nextIntraHopTime = determineHopTime(self.intraChainRate)
        elif self.nextInterHopTime == 0:
            self.nextInterHopTime = determineHopTime(self.interChainRate)
        else:
            print "self.nextIntraHopTime =", self.nextIntraHopTime
            print "self.nextInterHopTime =", self.nextInterHopTime
            raise SystemError('Neither hop time is zero - probably rounding error?')
        self.determineNextHopType()

        targetPosn = self.thiosInThisMolecule[hopTarget]
        self.position = targetPosn
        self.currentAtomIndex = hopTarget
        if plottingSubroutines == True:
            plotCarrier(self.position, self.thiosInThisMolecule, self.currentAtomIndex)
        

    def findBondedAtoms(self):
        bondedThiophenes = []

        for bond in self.thioThioBonds:
            if self.currentAtomIndex == bond[0]:
                bondedThiophenes.append(bond[1])
            elif self.currentAtomIndex == bond[1]:
                bondedThiophenes.append(bond[0])

        return bondedThiophenes
        


def plotCarrier(carrierLocation, moleculeBackbone, currentAtomIndex):
    moleculeToPlot = copy.deepcopy(moleculeBackbone)
    del moleculeToPlot[currentAtomIndex]
    x, y, z = zip(*moleculeToPlot.values())
    x = list(x)
    y = list(y)
    z = list(z)

    fig = P.figure()
    ax = p3.Axes3D(fig)
    for i in range(len(x)):
        ax.scatter(x[i], y[i], z[i], s = 20, c = 'b')

    ax.scatter(carrierLocation[0], carrierLocation[1], carrierLocation[2], s = 20, c = 'r')
    filelist = os.listdir('./')
    filenameCounter = 0
    for files in filelist:
        if ".png" in files:
            filenameCounter += 1
    filenameAddon = str(filenameCounter)
    while len(filenameAddon) < 3:
        filenameAddon = '0'+filenameAddon
    P.savefig('./test'+filenameAddon+'.png')



def run(bonds, trajectoryData, plotting=False):
    # Set a seed for reproducability within the random function calls
    R.seed(7960745560)
    # Decide whether to plot trajectory of the charge carrier (slow!)
    global plottingSubroutines
    plottingSubroutines = plotting
    # Sort the incoming trajectory into atom types by molecule number
    thioAtomsByMolecule, alk1AtomsByMolecule, alk2AtomsByMolecule = sortTrajByMolecule(trajectoryData)

    # Find a random thiophene to stat the charge carrier on and then initialise the carrier
    randomPos, molecule = randomPosition(thioAtomsByMolecule) # Also returns all of the thiophene atoms in this molecule
    hole = chargeCarrier(randomPos, thioAtomsByMolecule, molecule, bonds)
    for i in range(100):
        print "Executing hop", i
        hole.makeHop()

    
