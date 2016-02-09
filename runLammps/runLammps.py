import os
import numpy as np
import pylab as P
import csv
import subprocess

def getDatList(direc, mode):
    fileList = os.listdir(direc)
#    files = []
#    for element in fileList:
#        files.append(element[7:])
    datFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.dat'):
            if (mode == 'shrink'):
                if (mode in fileName) or (len(findIndex(fileName,'_')) == 5):
                    datFiles.append(str(fileName))
            elif (mode == 'Single'):
                if (len(findIndex(fileName,'_')) == 1):
                    datFiles.append(str(fileName))
            else:
                if (mode in fileName):
                    datFiles.append(str(fileName))
    
    # for folder in dataFolders:
    #     path = str(direc)+str(folder)+'/conv/'
    #     convList = os.listdir(path)
    #     for iterationFile in convList:
    #         dataFiles.append(str(path)+str(iterationFile))
    datFiles = sorted(datFiles)
    return datFiles


def selectFileToRun(validFiles):
    while True:
        print "\n---=== VALID FILES TO BE RUN ===---"
        if len(validFiles) == 0:
            print "ERROR: No valid '.dat' files found for this type!"
            return None
        for elementNo in range(len(validFiles)):
            print str(elementNo)+"):", validFiles[elementNo]
        print str(elementNo+1)+"): Exit runLammps.py"
        runThisFile = raw_input("Please select a file to run (integer, default = 0): ")
        if len(runThisFile) == 0:
            runThisFile = 0
        else:
            try:
                runThisFile = int(runThisFile)
            except:
                print "Please enter an integer between 0 and", len(validFiles)
                continue
        if (runThisFile < 0) or (runThisFile > len(validFiles)):
            print "Please enter an integer between 0 and", len(validFiles)
            continue
        elif runThisFile == len(validFiles):
            return None
        break
    return runThisFile



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
        inputTemplate[-2] = 'run              1000\n'
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
        inputFileLines[95] = 'dump           atomdump all custom 5 '+str(parameters[7])+str(dataFile)+'.lammpstrj id type mol x y z ix iy iz\n'
        inputFileLines[-2] = 'run               10000\n'
    elif parameters[6] == True: # Only one dumpstep at 16ps
        inputFileLines[95] = 'dump           atomdump all custom 100000 '+str(parameters[7])+str(dataFile)+'.lammpstrj id type mol x y z ix iy iz\n'
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
    return inputFileName, 1 #appendNumber


def configureSGE(lammpsInputFile, dataFile, appendNumber, parameters, procSlots, slurm=True):

    if slurm == True:
        SGETemplate = open('./templates/run.slurm', 'r')
    else:
        SGETemplate = open('./templates/run.sge', 'r')
    SGELines = SGETemplate.readlines()
    SGETemplate.close()

    # slashLoc = findIndex(lammpsInputFile, '/')
    # inputFileName = lammpsInputFile[slashLoc[-1]+1:]

#    folderLoc = lammpsInputFile[:slashLoc[-1]]+'/'
    
    GLoc = findIndex(SGELines[-1], '<')
    if slurm == True:
        SGELines[-1] = 'mpirun -np 8 ./lmp_ham < '+str(lammpsInputFile)[2:]+'\n'
    else:
        SGELines[-1] = 'mpirun -np $NSLOTS ./lmp_ham < '+str(lammpsInputFile)[2:]+'\n'
    
    if appendNumber != 1:
        if parameters[1] == True:
            if slurm == True:
                SGEName = './run_'+dataFile[:-4]+'_'+str(appendNumber)+'_shrink.slurm'
            else:
                SGEName = './run_'+dataFile[:-4]+'_'+str(appendNumber)+'_shrink.sge'
        elif parameters[2] == True:
            if slurm == True:
                SGEName = './run_'+dataFile[:-4]+'_'+str(appendNumber)+'_fastshrink.slurm'
            else:
                SGEName = './run_'+dataFile[:-4]+'_'+str(appendNumber)+'_fastshrink.sge'
        else:
            if slurm == True:
                SGEName = './run_'+dataFile[:-4]+'_'+str(appendNumber)+'.slurm'
            else:
                SGEName = './run_'+dataFile[:-4]+'_'+str(appendNumber)+'.sge'
        newFile = open(SGEName, 'w+')
        print "Writing SGE file as:", str(SGEName)
    else:
        if parameters[1] == True:
            if slurm == True:
                SGEName = './run_'+dataFile[:-4]+'_shrink.slurm'
            else:
                SGEName = './run_'+dataFile[:-4]+'_shrink.sge'
        elif parameters[2] == True:
            if slurm == True:
                SGEName = './run_'+dataFile[:-4]+'_fastshrink.slurm'
            else:
                SGEName = './run_'+dataFile[:-4]+'_fastshrink.sge'
        else:
            if slurm == True:
                SGEName = './run_'+dataFile[:-4]+'.slurm'
            else:
                SGEName = './run_'+dataFile[:-4]+'.sge'
        newFile = open(SGEName, 'w+')
        print "Writing SGE file as:", str(SGEName)

    newFile.writelines(SGELines)
    newFile.close()
    
    if procSlots > 1:
        if slurm == True:
            os.system("sbatch -p par6.q "+SGEName[2:])
        else:
            os.system("qsub -q par6.q -pe orte "+str(procSlots)+" "+SGEName[2:])
    else:
        if slurm == True:
            os.system("sbatch -p par6.q "+SGEName[2:])
        else:
            os.system("qsub -q seq6.q "+SGEName[2:])



def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
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
    direc = './'

    # Parameters of the form: [DoTestRun, DoSlowShrink, DoFastShrink, SimStartTemp, DoReduceTemperature, SimEndTemp, SingleChain, DumpSave, Timestep]
    parameters = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    simulationType = None


    while True:
        while True:
            print "\n---=== MODES THAT CAN BE EXECUTED BY RUNLAMMPS.PY ===---"
            print "0) Shrink a morphology"
            print "1) Heat a morphology"
            print "2) Cool a morphology"
            print "3) Equilibrate a morphology"
            print "4) Run a single chain"
            print "5) Exit program"
            mode = raw_input("Please select an execution mode (integer, default = 0): ")
            print "---====================================================---"
            if len(mode) == 0:
                mode = 0
                break
            else:
                try:
                    mode = int(mode)
                    break
                except:
                    print "Please enter an integer between 0 and 5."
                    continue
            if (mode < 0) or (mode > 5):
                print "Please enter an integer between 0 and 5."
                continue
            elif mode == 5:
                break
        if (mode == 0): # Fastshrink
            simulationType = 'Shrink'
            validDats = getDatList(direc, 'shrink')
            fileToRun = selectFileToRun(validDats)
            if fileToRun == None:
                break
            parameters[2] = True # Fastshrink
            parameters[3] = 290.0
        elif (mode == 1): # Heating
            simulationType = 'Heating'
            validDats = getDatList(direc, 'Heating')
            fileToRun = selectFileToRun(validDats)
            if fileToRun == None:
                break
            while True:
                temperature = raw_input("Please type in a temperature for the morphology to be run at (float, default 290.0): ")
                if len(temperature) == 0:
                    temperature = 290.0
                    break
                else:
                    try:
                        temperature = float(temperature)
                        break
                    except:
                        print "Please only use numerics."
                        continue
            parameters[4] = True
            parameters[3] = temperature
            parameters[5] = temperature
            # Rename the simulation to be run
            currentDatName = validDats[fileToRun][:-4]
            newDatName = currentDatName+'_'+str(int(temperature))+'K.dat'
            os.system("cp "+currentDatName+'.dat '+newDatName)
            validDats = [newDatName]
            fileToRun = 0
        elif (mode == 2): # Cooling
            simulationType = 'Cooling'
            validDats = getDatList(direc, 'Cooling')
            fileToRun = selectFileToRun(validDats)
            if fileToRun == None:
                break
            parameters[4] = True
            parameters[3] = float(validDats[fileToRun][-8:-5])
            parameters[5] = 290.0
        elif (mode == 3): # Equilibrate
            simulationType = 'Equilibrating'
            validDats = getDatList(direc, 'Equil')
            fileToRun = selectFileToRun(validDats)
            if fileToRun == None:
                break
            parameters[3] = 290.0
        elif (mode == 4): # Single chain
            simulationType = 'Single Chain'
            validDats = getDatList(direc, 'Single')
            fileToRun = selectFileToRun(validDats)
            if fileToRun == None:
                break
            parameters[3] = 290.0
            parameters[6] = True
        else:
            break

        parameters[0] = False
        parameters[7] = '/scratch/ghsk28/'
        parameters[8] = 4.0
        procSlots = 8

        print "\n---=== Current Simulation Parameters: ===---"
        print "Simulation Type =", simulationType
        print "Input Morphology =", validDats[fileToRun]
        print "Do small test run before submission =", bool(parameters[0])
        if (simulationType == 'Cooling') or (simulationType == 'Heating'):
            print "Reducing Temperature =", bool(parameters[4])
            print "Initial Temperature =", parameters[3]
            print "Final Temperature =", parameters[5]
        print "LAMMPSTRJ save location =", parameters[7]
        print "LAMMPS timestep =", parameters[8]
        print "Number of assigned processor slots =", procSlots
        print "---========================================---"
        while True:
            print "\n0) Accept current paramters and continue (default)"
            print "1) Modify parameters (not yet implemented)"
            paramOption = raw_input("Please choose whether to modify parameters: ")
            if (len(paramOption) == 0):
                paramOption = 0
                break
            else:
                try:
                    paramOption = int(paramOption)
                    break
                except:
                    paramOption = 0
                    break
            if (paramOption < 0) or (paramOption > 2):
                paramOption = 0
                break
            if paramOption == 0:
                print "Accepting these parameters..."
                break
            elif paramOption == 1:
                raise SystemError('NOT YET IMPLEMENTED!')
        

        lammpsInputFile, appendNumber = configureLammpsInput(validDats[fileToRun], parameters)
        if lammpsInputFile != 0:
            configureSGE(lammpsInputFile, validDats[fileToRun], appendNumber, parameters, procSlots)
        

    print "Battle control terminated."
