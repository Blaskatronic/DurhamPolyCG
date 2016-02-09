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
    morphologies = getEquilibratedMorphologies('./orcaOutputs')
    exitFlag = 0
    while exitFlag == 0:
        while True:
            print "\n---=== FILES THAT CAN BE RUN BY ANALYSEORCA.PY ===---"
            if len(morphologies) == 0:
                print "ERROR: No ORCA Outputs found in ./orcaOutputs. Please runOrca.py first!"
                exitFlag = 1
                break
            for elementNo in range(len(morphologies)):
                print str(elementNo)+"):", morphologies[elementNo]
            print str(elementNo+1)+"): Exit runAnalyseOrca.py"
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
                coresToUse = 0
                queue = None
                engine = None

            print "\n"
            t0 = T.time()
            print 'python analyseOrca.py -m '+str(morphologies[runThisFile])+' -c 0 -e None -q None'
            raise SystemError('STOP')
            subprocess.call('python analyseOrca.py -m '+str(morphologies[runThisFile])+' -c '+str(coresToUse)+' -e '+str(engine)+' -q '+str(queue), shell=True)
            # subprocess.call('python analyseOrca.py -m '+str(morphologies[runThisFile]), shell=True)
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
            print "analyseOrca.py calculations complete in %.1f %s." % (float(elapsedTime), str(timeunits))
            print "----------====================----------"

    print "Battle control terminated."
