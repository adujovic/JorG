#!/usr/bin/python

import sys
import argparse as ap

parser = ap.ArgumentParser(description='Convert from ascii text to \"textbytes\"')
parser.add_argument('--input', '-I', '-i', default=1,
                    help='input file')
parser.add_argument('--output', '-O', '-o', default=1,
                    help='output file')

flags = parser.parse_args()

if flags.__dict__['input'] == 1:
    text = [raw_input("Feed me input:\n\t")]
else:
    with open(flags.__dict__['input'],"r") as inFile:
        text = inFile.readlines()

output = ""
for arg in text:
    for c in arg:
        output+="%08d"%int(str(bin(ord(c)))[2:])

if flags.__dict__['output'] == 1:
    print('Output is:\n')
    print(output)
else:
    with open(flags.__dict__['output'],"w+") as outFile:
        outFile.write(output)
