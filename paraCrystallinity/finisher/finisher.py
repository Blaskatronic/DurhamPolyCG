import os
import numpy as np
import pylab as P
import random as R
import csv
import mpl_toolkits.mplot3d.axes3d as p3
import time as T
import copy
import itertools

def getFileList(direc):
    dataFiles = []
    fileList = os.listdir(direc)
#    files = []
#    for element in fileList:
#        files.append(element[7:])
    dataFiles = []
    for fileName in fileList:
        if (fileName[-9:-7] == '.o'):
            dataFiles.append(str(direc)+str(fileName))
    return dataFiles

def loadData(filename):
    timestepData = []
    paraData = []
    originalFile = open('./'+str(filename), 'r')
    dataLines = originalFile.readlines()
    originalFile.close()
    for line in dataLines:
        if "g = " in line:
            line = line[:-1]
            test = line.split(' ')
            paraData.append(float(test[-1]))
        if "dumpstep" in line:
            line = line[:-4]
            test = line.split(' ')
            timestepData.append(int(test[-1]))
        if "Traceback" in line:
            timestepData.pop(-1)
            break
    return timestepData, paraData


def plot(combinedTimestep, combinedData, combinedLabels, plottingOrder):
    colours = ['r','g','b','k','c','m','y']
    i = 0
    # Convert timesteps to realtime
    # The simulations are 1,000,000 timesteps with each timestep = 1fs
    # The dump files are every 2,000 timesteps so each dump represents 2ps
    # Paracrystallinity is taken every 50 dump steps = 0.1ns

    realTime = []
    fig = P.figure()
    ax = P.subplot(111)
    for dataNo in plottingOrder:
        # realTime = []
        # for j in range(len(combinedTimestep[dataNo])-1):
        #     realTime.append(0.1*j)
        # if len(combinedData[dataNo]) > len(realTime):
        #     combinedData[dataNo].pop(-1)
        ax.plot(combinedTimestep[dataNo], combinedData[dataNo], c=colours[i%len(colours)], label=combinedLabels[dataNo])
        print "---======---"
        print "Label =", combinedLabels[dataNo]
        print "Realtime =", combinedTimestep[dataNo]
        print "ParaData =", combinedData[dataNo]
        print "---======---"
#        P.plot(combinedTimestep[dataNo], combinedData[dataNo], c=colours[i%len(colours)], label=combinedLabels[dataNo])
        i += 1
    ax.set_xlabel('Simulation time (ns)')
    ax.set_ylabel('Paracrystallinity, g')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(flip(handles,4), flip(labels,4), loc = 'upper center', bbox_to_anchor=(0.5,1.1), ncol=4, fancybox = True)
    # bbox_to_anchor y-vals: Out of box (2 rows) = 1.145
    # 
    P.savefig('./combinedPara.png')



def findIndex(string, logical):
    '''This function returns the locations of an inputted character (logical) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == logical:
            locations.append(index)
        index += 1
    return locations


def flip(items, cols):
    return itertools.chain(*[items[i::cols] for i in range(cols)])
    


if __name__ == '__main__':
    combinedTimestep = []
    combinedData = []
    combinedLabels = []
    dataFiles = getFileList('./')
    print "Data Files", dataFiles
    i = 0
    plottingOrder = []
    for fileName in dataFiles:
        plottingOrder.append(0)

    dataFiles.sort()

#     for unfinishedFileName in dataFiles:
#         print "\n Examining:", unfinishedFileName

#         timestepData, paraData = loadData(unfinishedFileName)
#         fileName = unfinishedFileName[2:]
#         dotsList = findIndex(fileName, '.')
#         #name = fileName[:dotsList[0]]+'_'+fileName[dotsList[-1]+2:]
#         combinedTimestep.append(timestepData)
#         combinedData.append(paraData)
# #        combinedLabels.append(fileName[:dotsList[0]])
#         if ("continued" in fileName) or ("RTP" in fileName) or ("base" in fileName):
#             label = "290K"
#             plottingOrder[-1] = i
#         elif ("glass" in fileName) or ("Glass" in fileName):
#             label = "363K"
#             plottingOrder[-2] = i
#         elif ("anneal" in fileName) or ("Anneal" in fileName):
#             label = "423K"
#             plottingOrder[-3] = i
#         elif ("highT" in fileName) or ("HighT" in fileName):
#             label = "800K"
#             plottingOrder[-4] = i
#         combinedLabels.append(label)
#         i += 1

#     for unfinishedFileName in dataFiles:
#         print "\n Examining:", unfinishedFileName

#         timestepData, paraData = loadData(unfinishedFileName)
#         fileName = unfinishedFileName[2:]
#         dotsList = findIndex(fileName, '.')
#         #name = fileName[:dotsList[0]]+'_'+fileName[dotsList[-1]+2:]
#         combinedTimestep.append(timestepData)
#         combinedData.append(paraData)
# #        combinedLabels.append(fileName[:dotsList[0]])
#         if ("tol10" in fileName):
#             label = "10 A"
#             plottingOrder[-1] = i
#         elif ("tol09" in fileName):
#             label = "9 A"
#             plottingOrder[-2] = i
#         elif ("tol08" in fileName):
#             label = "8 A"
#             plottingOrder[-3] = i
#         elif ("tol07" in fileName):
#             label = "7 A"
#             plottingOrder[-4] = i
#         elif ("tol06" in fileName):
#             label = "6 A"
#             plottingOrder[-5] = i
#         elif ("tol05" in fileName):
#             label = "5 A"
#             plottingOrder[-6] = i
#         elif ("tol04" in fileName):
#             label = "4 A"
#             plottingOrder[-7] = i
#         elif ("tol03" in fileName):
#             label = "3 A"
#             plottingOrder[-8] = i
#         elif ("tol02" in fileName):
#             label = "2 A"
#             plottingOrder[-9] = i
#         elif ("tol01" in fileName):
#             label = "1 A"
#             plottingOrder[-10] = i
#         combinedLabels.append(label)
#         i += 1


#     for unfinishedFileName in dataFiles:
#         print "\n Examining:", unfinishedFileName

#         timestepData, paraData = loadData(unfinishedFileName)
#         fileName = unfinishedFileName[2:]
#         dotsList = findIndex(fileName, '.')
#         #name = fileName[:dotsList[0]]+'_'+fileName[dotsList[-1]+2:]
#         combinedTimestep.append(timestepData)
#         combinedData.append(paraData)
# #        combinedLabels.append(fileName[:dotsList[0]])
#         if ("test100" in fileName):
#             label = "100 Thios"
#             plottingOrder[-1] = i
#         elif ("tol" in fileName):
#             label = "All Thios"
#             plottingOrder[-2] = i
#         combinedLabels.append(label)
#         i += 1

    for unfinishedFileName in dataFiles:
        print "\n Examining:", unfinishedFileName

        timestepData, paraData = loadData(unfinishedFileName)
        fileName = unfinishedFileName[2:]
        dotsList = findIndex(fileName, '.')
        #name = fileName[:dotsList[0]]+'_'+fileName[dotsList[-1]+2:]
        combinedTimestep.append(timestepData)
        combinedData.append(paraData)
#        combinedLabels.append(fileName[:dotsList[0]])
        if ("dipTest" in fileName):
            label = "dipTest"
            plottingOrder[-1] = i
        combinedLabels.append(label)
        i += 1



    print "Plotting graph..."
    plot(combinedTimestep, combinedData, combinedLabels, plottingOrder)
    print "Graph plotted as ./combinedPara.png"
