#!/bin/bash

echo "starting tests:"
echo "python3 -m unittest tests.POSCARloader._unittest -v"
python3 -m unittest tests.POSCARloader._unittest -v
echo "                 ...done"

echo "python3 -m unittest tests.NearestNeighborsGenerator._unittest -v"
python3 -m unittest tests.NearestNeighborsGenerator._unittest -v
echo "                 ...done"

echo "python3 -m unittest tests.argv._unittest -v"
python3 -m unittest tests.argv._unittest -v
echo "                 ...done"
