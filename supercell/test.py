#!/usr/bin/python3

from aux.argv import options
import sys

a = options(*sys.argv)
print(a)
print(a('Wyckoffs'))
print(a('mask'))
print(a('neighbor'))
print(a('cutOff'))
