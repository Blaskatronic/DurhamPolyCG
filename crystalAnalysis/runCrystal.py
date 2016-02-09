import os
import numpy as np
import pylab as P
import csv
import subprocess


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



def findIndex(string, logical):
    '''This function returns the locations of an inputted character (logical) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == logical:
            locations.append(index)
        index += 1
    return locations


def configureSGE(lammpsInputFile, dataFile, queue):
    while (queue != 'par') and (queue != 'seq'):
        queue = raw_input("Queue name '"+str(queue)+"' invalid...please type 'par' or 'seq' to submit to valid Hamilton6 queue: ")
    SGETemplate = open('./templates/para.sge', 'r')
    SGELines = SGETemplate.readlines()
    SGETemplate.close()

    SGELines[-1] = SGELines[-1][:-1]

    SGELines[-1] += " -d "+str(dataFile)+" -l "+str(lammpsInputFile)+"\n"

    hyphenLocs = findIndex(dataFile, "_")
    dotLocs = findIndex(dataFile, ".")
    
    SGEName = "para_"+str(dataFile[:hyphenLocs[1]])+str(dataFile[hyphenLocs[-1]:dotLocs[-1]])+".sge"


    newFile = open(SGEName, 'w+')
    print "Writing SGE file as:", str(SGEName)

    newFile.writelines(SGELines)
    newFile.close()

    os.system("qsub -q "+str(queue)+"6.q "+SGEName)

if __name__ == "__main__":

    direc = './'
    validDats, validLammps = getFilesList(direc)
    usableFiles = zip(validDats, validLammps)
    usableFiles = sorted(usableFiles)

    exitFlag = 0

    dumpstep = -1
    hkl = [0, 1, 0]

    while exitFlag == 0:
        while True:
            print "\n---=== FILES THAT CAN BE RUN BY CRYSTAL.PY ===---"
            for elementNo in range(len(usableFiles)):
                print str(elementNo)+"):", usableFiles[elementNo]
            print str(elementNo+1)+"): Exit runCrystal.py"
            # print "Valid files =", zip(datFiles, lammpstrjFiles)
            runThisFile = raw_input("Please pick a file to run (integer, default = 0): ")
            if len(runThisFile) == 0:
                runThisFile = 0
            else:
                try:
                    runThisFile = int(runThisFile)
                except:
                    print "Please enter an integer between 0 and", len(usableFiles)
                    continue
            if (runThisFile < 0) or (runThisFile > len(usableFiles)):
                print "Please enter an integer between 0 and", len(usableFiles)
                continue
            elif runThisFile == len(usableFiles):
                print "Exiting Program..."
                exitFlag = 1
                break
            break
        if exitFlag == 0:
            print "\n---=== Current Paracrystallinity Calculation Settings: ===---"
            print "Datafile to use =", usableFiles[runThisFile][0]
            print "Corresponding LAMMPSTRJ file =", usableFiles[runThisFile][0]
            print "Dumpstep to be considered =", dumpstep
            print "---======================================================---"

            outputDir = './outputFiles/'
            name = str(usableFiles[runThisFile][1])[:-14]
            crystalPairsPickle = name+'_crystalPairs.pickle'
            fileList = os.listdir(outputDir)
            if crystalPairsPickle in fileList:
                print "Crystal pairs data found."
                pickleFound = True
                crystalCutOff = 0
            else:
                pickleFound = False
                print "No crystal pairs data found."
            # FIND THE CRYSTAL CUT-OFF IN THE .CSV FILE OR ASK THE USER AND THEN UPDATE THE CSV FILE
                crystalCutOff = float(raw_input("Please enter the crystal cut-off for this morphology: "))
                        
            # t0 = T.time()
            print 'python crystal.py -d '+str(usableFiles[runThisFile][0])+' -l '+str(usableFiles[runThisFile][1])+' -t '+str(dumpstep)+' -c '+str(crystalCutOff)+' -p '+str(int(pickleFound))
            subprocess.call('python crystal.py -d '+str(usableFiles[runThisFile][0])+' -l '+str(usableFiles[runThisFile][1])+' -t '+str(dumpstep)+' -c '+str(crystalCutOff)+' -p '+str(int(pickleFound)), shell=True)

    print "Battle control terminated."
