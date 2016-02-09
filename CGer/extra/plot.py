import os
import numpy as np
import pylab as P
import random as R
import csv
import mpl_toolkits.mplot3d.axes3d as p3
import time as T
import copy


if __name__ == '__main__':
    mean = 256
    yvals = []
    yvals2 = []
    for i in range(1000):
        yvals.append(R.gauss(256,30))

    for i in range(100,200,1):
        PDI = i/100.
        yvals2.append([])
        for j in range(1000000):
            yvals2[-1].append((R.random()**PDI)*PDI*mean)

    P.figure()
    P.hist(yvals2[-1], bins=100)
    P.show()
