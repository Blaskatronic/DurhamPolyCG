import multiprocessing as mp
import os
import time as T
import numpy as np
import itertools

def test(listObj, secondArg, thirdArg):
    T.sleep(2)
    print listObj, "second arg =", secondArg, "third arg =", thirdArg

def testStar(allArgs):
    return test(*allArgs)


if __name__ == '__main__':
    mp.freeze_support()
    listObj = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    count = 2
    pool = mp.Pool(processes=count)
    secondArg = 2
    thirdArg = 3
    pool.map(testStar, itertools.izip(listObj, itertools.repeat(secondArg), itertools.repeat(thirdArg)))
