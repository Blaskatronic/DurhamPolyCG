import random as R
import numpy as np
import time as T
import csv

if __name__ == '__main__':
    occupationMatrix = np.zeros((10,10,10))
    occupationMatrix2 = np.zeros((10,10,10))
    a = [5, 5, 5]
    t1 = T.time()
    for i in range(-2,3):
        for j in range(-2,3):
            for k in range(-2,3):
                occupationMatrix[(a[0]+i), (a[1]+j), (a[2]+k)] = 1
    t2 = T.time()
    occupationMatrix2[a[0]-2:a[0]+3, a[1]-2:a[1]+3, a[2]-2:a[2]+3] = 1
    t3 = T.time()

    print "First run =", t2-t1
    print "Second run =", t3-t2

    print occupationMatrix
    print "---======---"
    print occupationMatrix2
