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
import subprocess


def getEquilibratedMorphologies(direc):
    fileList = os.listdir(direc)
#    files = []
#    for element in fileList:
#        files.append(element[7:])
    equilibratedMorphologies = []
    for fileName in fileList:
        if ("Equil" in fileName):
            equilibratedMorphologies.append(str(fileName))
    return sorted(equilibratedMorphologies)


def checkDirectories(datFileName):
    dotsLoc = findIndex(datFileName, '.')
    actualName = datFileName[:dotsLoc[-1]]
    inputDir = os.listdir('./orcaInputs')
    outputDir = os.listdir('./orcaOutputs')
    if (actualName not in inputDir):
        if (actualName not in outputDir):
            print "Neither input nor output directories present!"
            os.makedirs('./orcaInputs/'+str(actualName))
            os.makedirs('./orcaOutputs/'+str(actualName))
            print "Correct directories created."
        else:
            print "Input directory not present!"
            os.makedirs('./orcaInputs/'+str(actualName))
            print "Correct directory created."
    else:
        print "Both input and output directories present and correct."


def configureSGE(morphology, cores, queue):
    SGEFilesNo = 1
    filesInDirectory = os.listdir('./')
    for fileName in filesInDirectory:
        if ".sge" == fileName[-4:]:
            SGEFilesNo += 1
    SGEFilesNo = str(SGEFilesNo)
    while len(SGEFilesNo) < 3:
        SGEFilesNo = '0'+SGEFilesNo
    SGEName = './runOrca'+SGEFilesNo+'.sge'
    templateSGE = open('./templates/serial.sge', 'r')
    templateLines = templateSGE.readlines()
    templateSGE.close()
    templateLines[-1] = 'python parallelRunOrca.py -m '+morphology+' -c '+cores+' -q '+queue

    newSGE = open(SGEName, 'w+')
    newSGE.writelines(templateLines)
    newSGE.close()

    return SGEName


def findFreeSlots(queue, engine):
    if engine == 'sge':
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
        return "<QUEUE NOT VISIBLE>"
    elif engine == 'slurm':
        slurmRequest = ['sinfo', '-o %10R%C']
        sinfoOutput = subprocess.Popen(slurmRequest, stdout=subprocess.PIPE).communicate()[0]
        sinfoOutput = sinfoOutput.split('\n')
        for lineNo in range(len(sinfoOutput)):
            if queue in sinfoOutput[lineNo]:
                queueLine = []
                tempLine = sinfoOutput[lineNo]
                tempLine = tempLine.split(' ')
                for element in tempLine:
                    if (len(element) != 0) and (element != "\n") and (queue not in element):
                        tempLine2 = element.split('/')
                        return int(tempLine2[1])
        return "<QUEUE NOT VISIBLE>"
    else:
        return "<ENGINE NOT FOUND>"



def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    return locations
   

def parallelSort(list1, list2):
    data = zip(list1, list2)
    data.sort()
    list1, list2 = map(lambda t: list(t), zip(*data))
    return list1, list2


if __name__ == '__main__':
    morphologies = getEquilibratedMorphologies('./orcaInputs')
    exitFlag = 0
    while exitFlag == 0:
        while True:
            print "\n---=== FILES THAT CAN BE RUN BY RUNORCAPARALLEL.PY ===---"
            if len(morphologies) == 0:
                print "ERROR: No ORCA Inputs found in ./orcaInputs. Please runSegments.py first!"
                exitFlag = 1
                break
            for elementNo in range(len(morphologies)):
                print str(elementNo)+"):", morphologies[elementNo]
            print str(elementNo+1)+"): Exit runOrca.py"
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

        if (exitFlag == 0):
            while True:
                submit = raw_input("Would you like to submit this job to Hamilton? (y/n, n runs it directly instead, y is default): ")
                if len(submit) == 0:
                    submit = 'y'
                try:
                    submit = str(submit)
                    if (submit != 'y') and (submit != 'n'):
                        raise Exception
                    else:
                        break
                except:
                    print "Please enter only y or n"

            if submit == 'y':
                parallel = True
                while True:
                    print "---=== AVAILABLE QUEUES ===---"
                    print "0) seq6.q (free slots on sge = "+str(findFreeSlots('seq6.q', 'sge'))+", free slots on slurm = "+str(findFreeSlots('seq6.q', 'slurm'))+")"
                    print "1) par6.q (free slots on sge = "+str(findFreeSlots('par6.q', 'sge'))+", free slots on slurm = "+str(findFreeSlots('par6.q', 'slurm'))+")"
                    print "2) seq5.q (free slots on sge = "+str(findFreeSlots('seq5.q', 'sge'))+", free slots on slurm = "+str(findFreeSlots('seq5.q', 'slurm'))+")"
                    print "3) par5.q (free slots on sge = "+str(findFreeSlots('par5.q', 'sge'))+", free slots on slurm = "+str(findFreeSlots('par5.q', 'slurm'))+")"
                    print "4) par4.q (free slots on sge = "+str(findFreeSlots('par4.q', 'sge'))+", free slots on slurm = "+str(findFreeSlots('par4.q', 'slurm'))+")"

                    try:
                        queue = raw_input("Please select a queue (integer, default = 0): ")
                        if len(str(queue)) == 0:
                            queue = 'seq6.q'
                            break
                        else:
                            if int(queue) == 0:
                                queue = 'seq6.q'
                                break
                            elif int(queue) == 1:
                                queue = 'par6.q'
                                break
                            elif int(queue) == 2:
                                queue = 'seq5.q'
                                break
                            elif int(queue) == 3:
                                queue = 'par5.q'
                                break
                            elif int(queue) == 4:
                                queue = 'par4.q'
                                break
                            else:
                                raise Exception
                    except:
                        raw_input("Invalid queue selection, please press return to try again...")

                while True:
                    print "---=== AVAILABLE SUBMISSION ENGINES ===---"
                    print "0) sge/current"
                    print "1) slurm/current"
                    try:
                        engine = raw_input("Please select an engine (integer, default = 0): ")
                        if len(str(engine)) == 0:
                            engine = 'sge'
                            break
                        else:
                            if int(engine) == 0:
                                engine = 'sge'
                                break
                            elif int(engine) == 1:
                                engine = 'slurm'
                                break
                            else:
                                raise Exception
                    except:
                        raw_input("Invalid engine selection, please press return to try again...")



                freeQueueSlots = findFreeSlots(queue, engine)

                while True:
                    coresToUse = raw_input("Please choose the number of cores to use for parallel ORCA resubmissions (default = 30): ")
                    if len(coresToUse) == 0:
                        coresToUse = 30
                    try:
                        coresToUse = int(coresToUse)
                    except:
                        print "Please type an integer that is greater than zero."
                        continue
                    if (coresToUse < 0):
                        print "Please type an integer that is greater than zero."
                        continue
                    if (coresToUse > freeQueueSlots):
                        print "Requested number of cores ("+str(coresToUse)+") is greater than the number of available cores on queue", queue, "("+str(freeQueueSlots)+")."
                        print "Jobs will likely sit in queue for a long time."
                        raw_input("Press return to continue...")
                    break
            else: #submit == 'n'
                parallel = False
                coresToUse = 0
                engine = None
                queue = None







        # if exitFlag == 0:
        #     while True:
        #         print "---=== AVAILABLE QUEUES ===---"
        #         print "0) seq6.q"
        #         print "1) par6.q"
        #         try:
        #             queue = raw_input("Please select a queue (integer, default = 0): ")
        #             if len(str(queue)) == 0:
        #                 queue = 'seq6.q'
        #                 break
        #             else:
        #                 if int(queue) == 0:
        #                     queue = 'seq6.q'
        #                     break
        #                 elif int(queue) == 1:
        #                     queue = 'par6.q'
        #                     break
        #                 else:
        #                     raise Exception
        #         except:
        #             print "Invalid queue selection, please try again."
        #             raw_input("Hit return to continue...")
        #     freeQueueSlots = findFreeSlots(queue)
        #     print "There are", freeQueueSlots, "free cores available in the", queue
        #     while True:
        #         if freeQueueSlots == 0:
        #             print "No slots available, running sequentially on current node using one core..."
        #             coresToUse = 0
        #             break
        #         coresToUse = raw_input("Please choose the number of cores to use for parallel ORCA submissions: ")
        #         try:
        #             coresToUse = int(coresToUse)
        #         except:
        #             print "Please type an integer between 0 and", freeQueueSlots
        #             continue
        #         if (coresToUse < 0):
        #             print "Please type an integer between 0 and", freeQueueSlots
        #             continue
        #         if (coresToUse > freeQueueSlots):
        #             print "Entered integer is greater than the number of free slots on seq6.q ("+str(freeQueueSlots)+")."
        #             print "Using", freeQueueSlots, "slots."
        #             coresToUse = freeQueueSlots
        #         break

            while True:
                t0 = T.time()
                print 'python parallelRunOrca.py -m '+str(morphologies[runThisFile])+' -c 0 -e None -q None'
                raise SystemError('STOP')
                subprocess.call('python parallelRunOrca.py -m '+str(morphologies[runThisFile])+' -c '+str(coresToUse)+' -e '+str(engine)+' -q '+str(queue), shell=True)
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
                print "ParallelRunOrca.py calculations complete in %.1f %s." % (float(elapsedTime), str(timeunits))
                print "----------====================----------"
                break

    print "Battle control terminated."
