#!/usr/bin/python

import sys
import argparse as ap

parser = ap.ArgumentParser(description='Convert from \"textbytes\" to ascii characters')
parser.add_argument('--input', '-I', '-i', default=1,
                    help='input file')
parser.add_argument('--output', '-O', '-o', default=1,
                    help='output file')

flags = parser.parse_args()

if flags.__dict__['input'] == 1:
    textbytes = [raw_input("Feed me input:\n\t")]
else:
    with open(flags.__dict__['input'],"r") as inFile:
        textbytes = inFile.readlines()

output = ""
for arg in textbytes:
    for byte in range(0,len(arg),8):
        char = 0
        for b in range(8):
            char += int(arg[byte+b])*(2**(7-b))
        output += "%c"%char

if flags.__dict__['output'] == 1:
    print('Output is:\n')
    print(output)
else:
    with open(flags.__dict__['output'],"w+") as outFile:
        outFile.write(output)

