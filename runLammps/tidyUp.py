import os
import numpy as np
import pylab as P
import random as R
import csv
import mpl_toolkits.mplot3d.axes3d as p3
import time as T
import copy



def tidyUpSingleChains(direc):
    fileList = os.listdir(direc)
    completedSims = []
    for element in fileList:
        if element[-14:] == '.dat.lammpstrj':
            dots = findIndex(element, '.')
            completedSims.append(element[:dots[0]])

    if len(completedSims) == 0:
        print "All single-chain simulations tidied! Job done."
        return

    for sim in completedSims:
#        os.system("ls *"+str(sim)+"*")
        os.system("mv *"+str(sim)+"* ./completed/singleChains/")

    print "All simulations tidied! Job done."
    return

def tidyUpOtherSims(direc):
    fileList = os.listdir(direc)
    completedSims = []
    for element in fileList:
        if element[-10] == '.lammpstrj':
            dots = findIndex(element, '.')
            completedSims.append(element[:dots[0]])

    if len(completedSims) == 0:
        print "All other simulations tidied! Job done."

    for sim in completedSims:
        os.system("mv *"+str(sim)+"* ./completed/miscSims/")
        
    
def findIndex(string, logical):
    '''This function returns the locations of an inputted character (logical) in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == logical:
            locations.append(index)
        index += 1
    return locations


if __name__ == '__main__':
    tidyUpSingleChains('./')
#    tidyUpOtherSims('./')
