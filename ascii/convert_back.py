#!/usr/bin/python

import argparse as ap

parser = ap.ArgumentParser(description='Convert from \"textbytes\" to ascii characters')
parser.add_argument('--input', '-I', '-i', default=None,
                    help='input file')
parser.add_argument('--output', '-O', '-o', default=None,
                    help='output file')
parser.add_argument('--back', '-B', action='store_true',
                    help='translation back')

flags = parser.parse_args()

try:
    with open(flags.__dict__['input'],"r") as inFile:
        TEXT = inFile.read()
except IOError:
    TEXT = raw_input("Feed me input:\n\t")

output = ""
if flags.__dict__['back']:
    for byte in range(0,len(TEXT),8):
        char = 0
        for b in range(8):
            char += int(TEXT[byte+b])*(2**(7-b))
        output += "%c"%char
else:
    for c in TEXT:
        output+="%08d"%int(str(bin(ord(c)))[2:])

try:
    with open(flags.__dict__['output'],"w+") as outFile:
        outFile.write(output)
except IOError:
    print('Output is:\n')
    print(output)
