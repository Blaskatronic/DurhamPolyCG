import multiprocessing as mp
import os
import time as T
import numpy as np
import itertools

def test(listObj, secondArg, thirdArg):
    print listObj, "second arg =", secondArg, "third arg =", thirdArg

def testStar(allArgs):
    return test(*allArgs)


if __name__ == '__main__':
    count = mp.cpu_count()
    print "I have found", count, "CPUs."
