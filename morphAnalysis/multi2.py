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
    mp.freeze_support()
    listObj = [1, 2, 3]
    count = 2
    for inputVar in listObj:
        newProc = mp.Process(target=test, args=(inputVar, 'Custard', 'cylons'))
        print newProc, newProc.is_alive()
        newProc.start()
        print newProc, newProc.is_alive()
        T.sleep(2)
        print newProc, newProc.is_alive()
        newProc.join()
        print newProc, newProc.is_alive()
        break
