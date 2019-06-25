#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../')
import time
from crun.crun import Crun

def main(**args):
    pass

if __name__ == '__main__':
    totalTracker  = -(time.time())
    tracker  = -(time.time())
    builder = Crun("hello/hello.cpp","hello/reference.cpp")
    builder("hello")
    builder("secondhello")
    del builder
    tracker += time.time()

    tracker  = -(time.time())
    builder = Crun("loop/loop.cpp")
    builder("loop",10)
    del builder
    tracker += time.time()

    tracker  = -(time.time())
    builder = Crun("sum/sum.cpp")
    a = builder("sum",11,12)
    print(a)
    del builder

    tracker += time.time()
    print("\t\tRuntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))

    totalTracker += time.time()
    print("\nTotal runtime: %02d:%02d:%02d"%(int(totalTracker/3600),int(totalTracker/60),int(totalTracker)))
