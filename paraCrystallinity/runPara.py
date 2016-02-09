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
            print "\n---=== FILES THAT CAN BE RUN BY PARA.PY ===---"
            for elementNo in range(len(usableFiles)):
                print str(elementNo)+"):", usableFiles[elementNo]
            print str(elementNo+1)+"): Exit runPara.py"
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
            while True:
                print "\n---=== Current Paracrystallinity Calculation Settings: ===---"
                print "Datafile to use =", usableFiles[runThisFile][0]
                print "Corresponding LAMMPSTRJ file =", usableFiles[runThisFile][0]
                print "Dumpstep to be considered =", dumpstep
                print "HKL axis to calculate paracrystallinity along =", hkl
                print "---======================================================---"

                while True:
                    print "\n0) Accept current parameters and continue (default)"
                    print "1) Change Dumpstep to be treated"
                    print "2) Modify HKL axis"
                    paramOption = raw_input("Please choose whether to modify parameters (integer): ")
                    if (len(paramOption) == 0):
                        paramOption = 0
                        break
                    else:
                        try:
                            paramOption = int(paramOption)
                        except:
                            paramOption = 0
                            break
                    if (paramOption < 0) or (paramOption > 2):
                        paramOption = 0
                        break
                    else:
                        break
                if paramOption == 0:
                    print "Accepting these parameters..."
                    break
                elif paramOption == 1:
                    newDumpstep = raw_input("Please key in new dumpstep number to consider (integer, default = -1): ")
                    if (len(newDumpstep) == 0):
                        newDumpstep = -1
                    else:
                        try:
                            newDumpstep = int(newDumpstep)
                        except:
                            print "New dumpstep not understood."
                            newDumpstep = -1
                    print "Dumpstep updated to:", newDumpstep
                    dumpstep = newDumpstep
                    continue
                elif paramOption == 2:
                    print "\n---=== HKL axis options ===---"
                    print "0) Pi-stacking direction (hkl = [0, 1, 0], default)"
                    print "1) Along-chain direction (hkl = [0, 0, 1])"
                    print "2) Alkyl-stacking direction (hkl = [1, 0, 0])"
                    print "---========================---"
                    hklOption = raw_input("Please choose an hkl axis (integer): ")
                    if (len(hklOption) == 0):
                        hklOption = 0
                    else:
                        try:
                            hklOption = int(hklOption)
                        except:
                            hklOption = 0
                    if (hklOption < 0) or (hklOption > 2):
                        hklOption = 0
                    

                    print "HKL axis updated to:",
                    if hklOption == 0:
                        print "Pi-stacking (hkl = [0, 1, 0])"
                        hkl = [0, 1, 0]
                    elif hklOption == 1:
                        print "Along-chain (hkl = [0, 0, 1])"
                        hkl = [0, 0, 1]
                    else:
                        print "Alkyl-stacking (hkl = [1, 0, 0])"
                        hkl = [1, 0, 0]
                    continue
                        
            # t0 = T.time()
            print 'python para.py -d '+str(usableFiles[runThisFile][0])+' -l '+str(usableFiles[runThisFile][1])+' -t '+str(dumpstep)+ ' -hkl "'+str(hkl)+'"'
            # raise SystemError('NO!')
            subprocess.call('python para.py -d '+str(usableFiles[runThisFile][0])+' -l '+str(usableFiles[runThisFile][1])+' -t '+str(dumpstep)+ ' -hkl "'+str(hkl)+'"', shell=True)
            # t1 = T.time()
            # elapsedTime = float(t1) - float(t0)
            # if elapsedTime < 60:
            #     timeunits = 'seconds.'
            # elif elapsedTime < 3600:
            #     elapsedTime /= 60.0
            #     timeunits = 'minutes.'
            # elif elapsedTime < 86400:
            #     elapsedTime /= 3600.0
            #     timeunits = 'hours.'
            # else:
            #     elapsedTime /= 86400.0
            #     timeunits = 'days.'
            # print "----------====================----------"
            # print "Segments.py calculations complete in %.1f %s." % (float(elapsedTime), str(timeunits))
            # print "----------====================----------"

    print "Battle control terminated."
