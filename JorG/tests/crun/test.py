#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../')
import time
from JorGpi.crun import Crun

def print_time(tracker):
    print("\t\tRuntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))

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
    print_time(tracker)

    tracker  = -(time.time())
    builder = Crun("loop/loop.cpp")
    builder("loop",10)
    del builder
    tracker += time.time()
    print_time(tracker)

    tracker  = -(time.time())
    builder = Crun("sum/sum.cpp")
    a = builder("sum",11,12)
    print(a)
    del builder
    tracker += time.time()
    print_time(tracker)

    options = {'extra_compile_args': ['-std=c++17','-O3','-Wall','-Wextra','-pedantic','-fopenmp'],
               'extra_link_args'   : ['-std=c++17','-fopenmp','-lm']}
    tracker  = -(time.time())
    builder = Crun("pi/pi.cpp",**options)
    builder("mandelbrot",78,44);
    del builder
    tracker += time.time()
    print_time(tracker)

    totalTracker += time.time()
    print_time(totalTracker)
