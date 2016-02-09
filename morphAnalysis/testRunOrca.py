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

if __name__ == '__main__':
    while True:
        T.sleep(1)
        output = subprocess.Popen("qstat", stdout=subprocess.PIPE).communicate()[0]
        print output
