import os
import numpy as np
import pylab as P
import csv
import subprocess

def getDatList(direc):
    fileList = os.listdir(direc)
#    files = []
#    for element in fileList:
#        files.append(element[7:])
    datFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.dat'):
            datFiles.append(str(fileName))
    
    # for folder in dataFolders:
    #     path = str(direc)+str(folder)+'/conv/'
    #     convList = os.listdir(path)
    #     for iterationFile in convList:
    #         dataFiles.append(str(path)+str(iterationFile))
    # dataFiles = sorted(dataFiles)
    return datFiles

def configureLammpsInput(dataFile, parameters):
    print "Configuring LAMMPS input file for", dataFile
    inputHandle = open('./templates/run.in', 'r')
    inputTemplate = inputHandle.readlines()
    inputHandle.close()

    dataHandle = open('./'+str(dataFile), 'r')
    dataLines = dataHandle.readlines()
    dataHandle.close()

    description = dataLines[0]
    inputTemplate[0] = '# '+description

    inputTemplate[17] = 'read_data       '+str(dataFile)+'\n'     # Read_data line
    inputTemplate[18] = '#read_restart    '+str(dataFile)+'\n'    # Read_restart line

#    temperature = str(float(raw_input('Please input the required temperature of the simulation: ')))
    startTemperature = str(parameters[3])
    if parameters[4] == True: # doChangeTemperature
        endTemperature = str(parameters[5])
    else:
        endTemperature = str(parameters[3])

    # print "Temp:", inputTemplate[-18] # Temp Line
    # print "Compress:", inputTemplate[-17] # Compress line
    # print "Vel:", inputTemplate[-12] # Velocity Line (also temp)
    # print "First:", inputTemplate[-10] # First dump line
    # print "Second:", inputTemplate[-9] # Second dump line
    # print "Restart:", inputTemplate[-7] # Restart line
    # print "Run:", inputTemplate[-2] # Run line

    tempLine = inputTemplate[-18]
    test = tempLine.split(' ')
    newTempLine = ''
    for elementNo in range(len(test)):
        if (test[elementNo] == 'dummytemp'):
            test[elementNo] = startTemperature
        elif (test[elementNo] == 'dummytemp2'):
            test[elementNo] = endTemperature
        newTempLine += test[elementNo]+' '
    inputTemplate[-18] = newTempLine


    velLine = inputTemplate[-12]
    test = velLine.split(' ')
    newVelLine = ''
    for elementNo in range(len(test)):
        if (test[elementNo] == 'dummytemp'):
            test[elementNo] = startTemperature
        # elif (test[elementNo] == 'dummytemp2'): # Not actually necessary - it's only the initial velocity
        #     test[elementNo] = endTemperature
        newVelLine += test[elementNo]+' '
    inputTemplate[-12] = newVelLine

    firstDump = inputTemplate[-10]
    test = firstDump.split(' ')
    newFirstLine = ''
    for elementNo in range(len(test)):
        if (test[elementNo] == 'dummy.lammpstrj'):
            if parameters[1] == True: # Do shrink
                test[elementNo] = str(parameters[7])+str(dataFile)+'_shrink.lammpstrj'
            elif parameters[2] == True: # Do fast shrink
                test[elementNo] = str(parameters[7])+str(dataFile)+'_fastshrink.lammpstrj'
            else:
                test[elementNo] = str(parameters[7])+str(dataFile)+'.lammpstrj'
        newFirstLine += test[elementNo]+' '
    inputTemplate[-10] = newFirstLine

    secondDump = inputTemplate[-9]
    test = secondDump.split(' ')
    newSecondLine = ''
    for elementNo in range(len(test)):
        if (test[elementNo] == 'dummy.lammpstrj'):
            if parameters[1] == True: # Do shrink
                test[elementNo] = str(parameters[7])+str(dataFile)+'_shrink.lammpstrj'
            elif parameters[2] == True: # Do fast shrink
                test[elementNo] = str(parameters[7])+str(dataFile)+'_fastshrink.lammpstrj'
            else:
                test[elementNo] = str(parameters[7])+str(dataFile)+'.lammpstrj'
        newSecondLine += test[elementNo]+' '
    inputTemplate[-9] = newSecondLine

    restart = inputTemplate[-7]
    test = restart.split(' ')
    newRestartLine = ''
    for elementNo in range(len(test)):
        if (test[elementNo] == 'dummy.restart'):
            test[elementNo] = str(dataFile)+'.restart'
        newRestartLine += test[elementNo]+' '
    inputTemplate[-7] = newRestartLine

    # Set timestep:
    inputTemplate[100] = inputTemplate[100][:-6]+str(float(parameters[8]))+'  \n'


#### UNCOMMENT TO DO TEST RUNS
#     # Configure the run.in to a small test run of 100 timesteps for the check

    if parameters[0] == True: # Do test run
        inputTemplate[-2] = 'run              100\n'
        inputTemplate[-10] = '#'+inputTemplate[-10]
        inputTemplate[-9] = inputTemplate[-9][1:]

    underscoreLoc = findIndex(dataFile, '_')
# #    folderName = dataFile[:underscoreLoc[0]]

#  #   print "Folder Name =", folderName
    
    # appendNumber = checkFilesWithSameName(dataFile, '.in')
    # if appendNumber != 1:
    #     inputFileName = './'+dataFile[:-4]+'_'+str(appendNumber)+'.in'
    #     newFile = open(inputFileName, 'w+')
    #     print "Writing LAMMPS input file as:", str(inputFileName)
    # else:
    inputFileName = './'+dataFile[:-4]+'.in'
    if parameters[1] == True: # Do shrink
        inputFileName = './'+dataFile[:-4]+'_shrink.in'
    if parameters[2] == True: # Do fast shrink
        inputFileName = './'+dataFile[:-4]+'_fastshrink.in'

    newFile = open(inputFileName, 'w+')
    print "Writing LAMMPS input file as:", str(inputFileName)

    newFile.writelines(inputTemplate)
    newFile.close()

    if parameters[0] == True:
        testComplete = runLammpsTest(inputFileName)
        if testComplete == 0:
            print "Check input data file for molecule overlap, errors, bugs etc."
            return 0, 0
        os.system("rm "+str(parameters[7])+str(inputFileName[:-3])+".lammpstrj")
        os.system("rm ./log.lammps")
        print "Setting up main run..."
    
    inputFileRead = open(inputFileName, 'r')
    inputFileLines = inputFileRead.readlines()
    inputFileRead.close()
    if parameters[2] == True: # Do a fast shrink
        inputFileLines[-2] = 'run               100000\n'
    elif parameters[6] == True: # Only one dumpstep at 16ps
        inputFileLines[95] = 'dump           atomdump all custom 100000 '+str(dataFile)+'.lammpstrj id type mol x y z ix iy iz\n'
        inputFileLines[-2] = 'run              1000000\n'
    else:
        inputFileLines[-2] = 'run              1000000\n'
    if parameters[0] == True: # Did a test run, so restore the dumps
        inputFileLines[-10] = inputFileLines[-10][1:]
        inputFileLines[-9] = '#'+inputFileLines[-9]
    if (parameters[1] == True) or (parameters[2] == True): # Shrink volume (so uncomment)
        tempLine = inputFileLines[-17]
        newTempLine = tempLine[2:]
        inputFileLines[-17] = newTempLine
    inputFileWrite = open(inputFileName, 'w')
    inputFileWrite.writelines(inputFileLines)
    inputFileWrite.close()
    return inputFileName, 1#appendNumber

def configureSGE(lammpsInputFile, dataFile, appendNumber, parameters, procSlots):
    SGETemplate = open('./templates/run.sge', 'r')
    SGELines = SGETemplate.readlines()
    SGETemplate.close()

    # slashLoc = findIndex(lammpsInputFile, '/')
    # inputFileName = lammpsInputFile[slashLoc[-1]+1:]

#    folderLoc = lammpsInputFile[:slashLoc[-1]]+'/'
    
    GLoc = findIndex(SGELines[-1], '<')
    SGELines[-1] = 'mpirun -np $NSLOTS ./lmp_ham < '+str(lammpsInputFile)[2:]+'\n'
    
    if appendNumber != 1:
        if parameters[1] == True:
            SGEName = './run_'+dataFile[:-4]+'_'+str(appendNumber)+'_shrink.sge'
        elif parameters[2] == True:
            SGEName = './run_'+dataFile[:-4]+'_'+str(appendNumber)+'_fastshrink.sge'
        else:
            SGEName = './run_'+dataFile[:-4]+'_'+str(appendNumber)+'.sge'
        newFile = open(SGEName, 'w+')
        print "Writing SGE file as:", str(SGEName)
    else:
        if parameters[1] == True:
            SGEName = './run_'+dataFile[:-4]+'_shrink.sge'
        elif parameters[2] == True:
            SGEName = './run_'+dataFile[:-4]+'_fastshrink.sge'
        else:
            SGEName = './run_'+dataFile[:-4]+'.sge'
        newFile = open(SGEName, 'w+')
        print "Writing SGE file as:", str(SGEName)

    newFile.writelines(SGELines)
    newFile.close()

    if procSlots > 1:
        os.system("qsub -q par6.q -pe orte "+str(procSlots)+" "+SGEName[2:])
    else:
        os.system("qsub -q seq6.q "+SGEName[2:])



def findIndex(string, logical):
    '''This function returns the locations of an inputted character (logical) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == logical:
            locations.append(index)
        index += 1
    return locations
    

def checkFilesWithSameName(fileName, fileExt):
    fileList = os.listdir('./')
    dotsList = findIndex(fileName, '.')
    noOfFilesWithSameName = 1
    for element in fileList:
        if ((str(dataFile[:-4])+str(fileExt) in element) or (str(dataFile[:-4])+'_' in element)) and (fileExt in element):
            noOfFilesWithSameName += 1
    return noOfFilesWithSameName


def runLammpsTest(inputFileName):
#    print "./"+folderName+"/lmp_ham < "+inputFileName
    os.system("./templates/lmp_ham < "+inputFileName)
    while True:
        print "\nDid the simulation work? Type 'y' or 'n':"
        testComplete = str(raw_input(' '))
        if testComplete == 'y':
            return 1
        elif testComplete == 'n':
#            os.system("rm "+inputFileName)
            return 0
    



if __name__ == "__main__":
    # Copy DAT into this polymers folder:

    # Get the DAT and put it in the right place
    # Get the temperature required
    # Change the run.in and run.sge files to correspond to the DAT
    # Run a small LAMMPS run
    # Submit full LAMMPS to Hamilton if the small test works

    validDats = getDatList('./')
    validDats = sorted(validDats)
    print "Valid .dat files:", validDats


    doTestRun = True
    singleChain = False
#    doShrinkVolume = False
    doFastShrink = False
    doReduceTemperature = False

    simulationStartTemperature = 290.0 # RTP: 290, P3HT glass transition: 363, Common annealing at: 423, "HighT": 800
    simulationEndTemperature = 290.0 # RTP: 290, P3HT glass transition: 363, Common annealing at: 423, "HighT": 800
    
    timestep = 4.0
    procSlots = 8


    dumpSave = '/scratch/ghsk28/'
    parameters = [doTestRun, 0, doFastShrink, simulationStartTemperature, doReduceTemperature, simulationEndTemperature, singleChain, dumpSave, timestep]

    foundFile = 0


###############################################################################
# CUSTOM SIMULATING
###############################################################################
    datPart = 'Mw_70_PDI_1.6_con_40_Heating_290K.dat'



    # Find and submit ALL files with datPart in name:
    filesToRun = [] # Of the form: [[name, simStartTemp, simEndTemp], ...]

    while True:
        for dataFile in validDats:
            if datPart in dataFile:
                if doReduceTemperature == True:
                    if "Cooling" in dataFile:                    
                        foundFile = 1
                        parameters[3] = float(dataFile[-8:-5]) # simulationStartTemperature
                        parameters[5] = 290.0                  # simulationEndTemperature
                        filesToRun.append([dataFile, parameters[3], parameters[5]])
                    else:
                        raise SystemError('doReduceTemperature == True, but "Cooling" not in dataFile name')
                else:
                    if "Heating" in dataFile:
                        foundFile = 1
                        parameters[3] = float(dataFile[-8:-5]) # simulationStartTemperature
                        parameters[5] = float(dataFile[-8:-5]) # simulationEndTemperature
                        filesToRun.append([dataFile, parameters[3], parameters[5]])
                    elif "Equil" in dataFile:
                        foundFile = 1
                        parameters[3] = float(290.0)
                        parameters[5] = float(290.0)
                        filesToRun.append([dataFile, parameters[3], parameters[5]])
                    else:
                        raise SystemError('doReduceTemperature == False, but "Heating" or "Equil" not in dataFile name')

        if foundFile == 1:
            print "Using datafiles:", filesToRun
            break
        else:
            raise SystemError('No datafile found with "'+str(datPart)+'" in filename')

    for fileStats in filesToRun:
        dataFile = fileStats[0]
        parameters[3] = fileStats[1]
        parameters[5] = fileStats[2]
        lammpsInputFile, appendNumber = configureLammpsInput(dataFile, parameters)

        if lammpsInputFile != 0:
            configureSGE(lammpsInputFile, dataFile, appendNumber, parameters, procSlots)


    print "All done!"
    
###############################################################################




