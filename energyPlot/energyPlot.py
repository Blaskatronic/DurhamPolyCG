import os
import numpy as np
import pylab as P
import csv

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

def readData(dataFile):
    print "Attempting to read in data..."

    dataHandle = open(str(dataFile))
    rawData = dataHandle.readlines()
    dataHandle.close()

    # remove Header and Footer
    rawData = rawData[23:-21]

#    print rawData

    ###########################################################################################################################
    ##                                                                                                                       ##
    ## rawData is now in the form:                                                                                           ##
    ## [[ Step, CPU, Temp, PotEng, TotEng, Press, Volume, E_pair, E_bond, E_angle, E_dihed, E_impro ], FOR EACH TIMESTEP ]    ##
    ##                                                                                                                       ##
    ###########################################################################################################################


    steps = []
    potEng = []
    totEng = []
    EPair = []
    EBond = []
    EAngle = []
    EDihed = []
    EImpro = []
    temp = []

    for line in rawData:
        test = line.split(' ')
        removes = []
        for position in range(len(test)):
            if len(test[position]) == 0:
                removes.append(position)
        removes = sorted(removes, reverse=True)
        for position in removes:
            test.pop(position)
        test.pop(-1)
        try:
            for valueNo in range(len(test)):
                test[valueNo] = float(test[valueNo])
        except:
            continue
        try:
            steps.append(test[0])
            temp.append(test[2])
            potEng.append(test[3])
            totEng.append(test[4])
            EPair.append(test[7])
            EBond.append(test[8])
            EAngle.append(test[9])
            EDihed.append(test[10])
            EImpro.append(test[11])
        except:
            print "Did the simulation fail? Breaking read..."
            break

    return steps, temp, potEng, totEng, EPair, EBond, EAngle, EDihed, EImpro


def plotGraph(fileName, steps, temp, potEng, totEng, EPair, EBond, EAngle, EDihed, EImpro):
    print "Plotting Total Energy Graph..."
    P.figure()
    P.plot(steps, potEng, c='red', label='PotEng')
    P.plot(steps, totEng, c='blue', label='TotEng')
    P.ylabel('Energy, kcal/mol')
    P.xlabel('Timestep')
    P.legend(('PotEng', 'TotEng'))
    P.xlim((0,10000))
    P.savefig('./outputGraphs/'+str(fileName)+'_TotEnergy.png')
    P.clf()
    print "Graph saved as ./outputGraphs/"+str(fileName)+"_TotEnergy.png"


    print "Plotting Combined Energy Graph..."
    P.figure()
    P.plot(steps, EPair, c ='red', label='Pair')
    P.plot(steps, EBond, c ='blue', label='Bond')
    P.plot(steps, EAngle, c ='green', label='Angle')
    P.plot(steps, EDihed, c = 'yellow', label='Dihed')
    P.plot(steps, EImpro, c = 'black', label='Impro')
    P.ylabel('Energy, kcal/mol')
    P.xlabel('Timestep')
    P.legend(('Pair', 'Bond', 'Angle', 'Dihed', 'Impro'))
    P.xlim((0,np.max(steps)))
    P.savefig('./outputGraphs/'+str(fileName)+'_CombEnergy.png')
    P.clf()
    print "Graph saved as ./outputGraphs/"+str(fileName)[2:]+"_CombEnergy.png"

    print "Plotting Temperature Decay Curve..."
    P.figure()
    P.plot(steps, temp, c='red')
    P.ylabel('Temperature, K')
    P.xlabel('Timestep')
    P.xlim((0, np.max(steps)))
    P.savefig('./outputGraphs/'+str(fileName)+'_Temp.png')
    P.clf()
    print "Graph saved as ./outputGraphs/"+str(fileName)[2:]+"_Temp.png"

def plotGraphCombined(combinedSteps, combinedTotalEnergy, combinedLabels, plottingOrder):
    colours = ['r','g','b','k','c','m','y']
    i = 0
    P.figure()
    for dataNo in plottingOrder:
#    for dataNo in range(len(combinedTotalEnergy)):
        P.plot(combinedSteps[dataNo], combinedTotalEnergy[dataNo], c=colours[i%len(colours)], label=combinedLabels[dataNo])
        i += 1
    P.xlabel('Timestep')
    P.ylabel('Energy, kcal/mol')
    P.legend(loc = 'upper center', bbox_to_anchor=(0.5,1.05), ncol=4, fancybox = True)
    P.savefig('./outputGraphs/combinedEnergy.png')


def findIndex(string, logical):
    '''This function returns the locations of an inputted character (logical) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == logical:
            locations.append(index)
        index += 1
    return locations
    


def makeCSV(fileName, steps, potEng):
    csvFile = csv.writer(open('./outputCSVs/'+str(fileName)+'_potEng.csv', 'w+'), delimiter = ',')
    for stepNo in range(len(steps)):
        csvFile.writerow([steps[stepNo], potEng[stepNo]])
    print "CSV file saved as ./outputCSVs/"+str(fileName)+"_potEng.csv"



if __name__ == "__main__":
    dataFiles = getFileList('./')
    print "Data Files", dataFiles
    combinedSteps = []
    combinedTotalEnergy = []
    combinedLabels = []
    i = 0
    plottingOrder = []
    for fileName in dataFiles:
        plottingOrder.append(0)

    for fileName in dataFiles:
         print "\nExamining:", fileName
         steps, temp, potEng, totEng, EPair, EBond, EAngle, EDihed, EImpro = readData(fileName)
         plotGraph(fileName, steps, temp, potEng, totEng, EPair, EBond, EAngle, EDihed, EImpro)
         makeCSV(fileName, steps, potEng)
#          combinedSteps.append(steps)
#          combinedTotalEnergy.append(totEng)
#          name = fileName[2:]
# #         dotsList = findIndex(name, '.')
# #         label = name[:dotsList[0]]
#          if "290K" in name:
#              label = "290K"
#              plottingOrder[8] = i
#          elif "363K" in name:
#              label = "363K"
#              plottingOrder[7] = i
#          elif "423K" in name:
#              label = "423K"
#              plottingOrder[6] = i
#          elif "473K" in name:
#              label = "473K"
#              plottingOrder[5] = i
#          elif "673K" in name:
#              label = "673K"
#              plottingOrder[4] = i
#          elif "1073K" in name:
#              label = "1073K"
#              plottingOrder[3] = i
#          elif "1273K" in name:
#              label = "1273K"
#              plottingOrder[2] = i
#          elif "2273K" in name:
#              label = "2273K"
#              plottingOrder[1] = i
#          elif "4273K" in name:
#              label = "4273K"
#              plottingOrder[0] = i
             
#          combinedLabels.append(label)
#          i += 1

#     print "\nPlotting combined graph..."
#     plotGraphCombined(combinedSteps, combinedTotalEnergy, combinedLabels, plottingOrder)
#     print "Graph saved as ./outputGraphs/combinedEnergy.png"

    print "All done!"


