import os
import numpy as np
import pylab as P
import time as T
import sys
import argparse
import subprocess

def getFilesList(direc):
    fileList = os.listdir(direc)
    morphologyFiles = []
    for fileName in fileList:
        if ("Mw" in fileName):
            morphologyFiles.append(str(fileName))
    # Now check that the corresponding dat files are present
    homeDir = os.listdir('./')
    popList = []
    for morphologyNo in range(len(morphologyFiles)):
        if str(morphologyFiles[morphologyNo])+'.dat' not in homeDir:
            popList.append(morphologyNo)
    popList.sort(reverse = True)
    for element in popList:
        morphologyFiles.pop(element)
    return sorted(morphologyFiles)


def findFreeSlots(queue):
    qfreeOutput = subprocess.Popen("qfree", stdout=subprocess.PIPE).communicate()[0]
    qfreeOutput = qfreeOutput.split('\n')
    for lineNo in range(len(qfreeOutput)):
        if queue in qfreeOutput[lineNo]:
            # This is the line we need
            queueLine = []
            tempLine = qfreeOutput[lineNo].split(' ')
            for element in tempLine:
                if (len(element) != 0) and (element != "\n"):
                    queueLine.append(element)
            # Output looks something like [QUEUE, PE, UsedFor/By, FREE, USED, TOTAL]
            # UsedFor/By might have multiple spaces and text in though depending on the queue
            # So instead of calling element [3], we should pick element [-3] to get the number of free slots
            return int(queueLine[-3])
    print "Correct queueline for '"+str(queue)+"' not found. Running on current node only..."
    return 0

def configureSGE(morphology, time, carriers, cores, queue, plotting):
    SGEFilesNo = 1
    filesInDirectory = os.listdir('./')
    for fileName in filesInDirectory:
        if fileName[-4:] == ".sge":
            SGEFilesNo += 1
    SGEFilesNo = str(SGEFilesNo)
    while len(SGEFilesNo) < 3:
        SGEFilesNo = '0'+SGEFilesNo
    SGEName = './runCT'+SGEFilesNo+'.sge'
    templateSGE = open('./MorphologyAnalysis/templates/serial.sge', 'r')
    templateLines = templateSGE.readlines()
    templateSGE.close()
    templateLines[-1] = 'python ./CTmod.py -m '+str(morphologies[runThisFile])+' -t '+str(simulationTime)+' -c '+str(carriers)#+' -n '+str(coresToUse)+' -q '+str(queue)+' -p '+str(plotting)

    newSGE = open(SGEName, 'w+')
    newSGE.writelines(templateLines)
    newSGE.close()

    return SGEName


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



def trimLammpstrj(filename):
    timestepData = []
    timestepNos = []
    simVolData = []
    print "Reading in trajectory file:", filename
    originalFile = open('./'+str(filename), 'r')
    trajData = originalFile.readlines()
    originalFile.close()

    newTimestepLines = []

    for lineNo in range(len(trajData)):
        if "ITEM: TIMESTEP" in trajData[lineNo]:
            newTimestepLines.append(lineNo)

    trimmedTrajData = trajData[newTimestepLines[-1]:]
    

    trimmedFile = open('./'+str(filename)+'.trim', 'w+')
    trimmedFile.writelines(trimmedTrajData)
    trimmedFile.close()


    return './'+str(filename)+'.trim'



if __name__ == "__main__":

    # Throw in the same stuff from the other programs to pick which morphology to run. For now just run the one thing that's in there
    simulationTimes = np.logspace(-14,-4,11)
    carriers = 100

    morphologies = getFilesList('./preparedMorphologies')
    exitFlag = 0
    while exitFlag == 0:
        while True:
            print "\n---=== FILES THAT CAN BE RUN BY CTMOD.PY ===---"
            if len(morphologies) == 0:
                print "ERROR: No Morphologies found in ./preparedMorphologies with matching .dat files for the morphology in this directory './'. Please copy the results from analyseOrca.py to the correct locations first!"
                exitFlag = 1
                break
            for elementNo in range(len(morphologies)):
                print str(elementNo)+"):", morphologies[elementNo]
            print str(elementNo+1)+"): Exit runCT.py"
            # print "Valid files =", zip(datFiles, lammpstrjFiles)
            runThisFile = raw_input("Please pick a file to run (integer, default = 0): ")
            if len(runThisFile) == 0:
                runThisFile = 0
            else:
                try:
                    runThisFile = int(runThisFile)
                except:
                    print "Please enter an integer between 0 and", len(morphologies)
                    continue
            if (runThisFile < 0) or (runThisFile > len(morphologies)):
                print "Please enter an integer between 0 and", len(morphologies)
                continue
            elif runThisFile == len(morphologies):
                print "Exiting Program..."
                exitFlag = 1
                break
            break
        if exitFlag == 0:
            while True:
                print "---=== AVAILABLE QUEUES ===---"
                print "0) seq6.q"
                print "1) par6.q"
                try:
                    queue = raw_input("Please select a queue (integer, default = 0): ")
                    if len(str(queue)) == 0:
                        queue = 'seq6.q'
                        break
                    else:
                        queue = int(queue)
                        if queue == 0:
                            queue = 'seq6.q'
                            break
                        elif queue == 1:
                            queue = 'par6.q'
                            break
                        else:
                            raise Exception
                except:
                    print "Invalid queue selection, please try again."
                    raw_input("Hit return to continue...")
            freeQueueSlots = findFreeSlots(queue)
            print "There are", freeQueueSlots, "free cores available in the", queue
            coresToUse = 1


            while True:
                submit = raw_input("Would you like to submit this job to Hamilton? (y/n, n runs it directly instead): ")
                if len(submit) == 0:
                    submit = 'n'
                try:
                    submit = str(submit)
                    if submit == 'y':
                        for simulationTime in simulationTimes:
                            SGEName = configureSGE(morphologies[runThisFile], simulationTime, carriers, coresToUse, queue, 0)
                            if queue == "seq6.q":
                                print "qsub -q "+queue+" "+SGEName[2:]
                                os.system("qsub -q "+queue+" "+SGEName[2:])
                            elif queue == "par6.q":
                                print "qsub -q "+queue+" "+SGEName[2:]
                                os.system("qsub -q "+queue+" -pe orte 1 "+SGEName[2:])
                            else:
                                print "Queue =", queue
                                print "Queue not correctly hardcoded. No submissions have occured."
                                break
                        break
                    elif submit == 'n':
                        t0 = T.time()
                        for simulationTime in simulationTimes:
                            subprocess.call('python CTmod.py -m '+str(morphologies[runThisFile])+' -t '+str(simulationTime)+' -c '+str(carriers)+' -n '+str(coresToUse)+' -q '+str(queue)+' -p 0', shell=True)
                        t1 = T.time()
                        elapsedTime = float(t1) - float(t0)
                        if elapsedTime < 60:
                            timeunits = 'seconds.'
                        elif elapsedTime < 3600:
                            elapsedTime /= 60.0
                            timeunits = 'minutes.'
                        elif elapsedTime < 86400:
                            elapsedTime /= 3600.0
                            timeunits = 'hours.'
                        else:
                            elapsedTime /= 86400.0
                            timeunits = 'days.'
                        print "----------====================----------"
                        print "CTmod.py calculations complete in %.1f %s." % (float(elapsedTime), str(timeunits))
                        print "----------====================----------"
                        break
                    else:
                        raise Exception
                except:
                    print "Please enter only y or n."
                    raw_input("Press return to continue...")








    print "Battle control terminated."






        
    
