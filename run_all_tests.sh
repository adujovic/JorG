#!/bin/bash
cd JorG

ERR=0

echo "starting tests:"
echo "python3 -m unittest tests.POSCARloader._unittest -v"
python3 -m unittest tests.POSCARloader._unittest -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi
echo "                 ...done"

echo "python3 -m unittest tests.NearestNeighborsGenerator._unittest -v"
python3 -m unittest tests.NearestNeighborsGenerator._unittest -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi
echo "                 ...done"

echo "python3 -m unittest tests.argv._unittest -v"
python3 -m unittest tests.argv._unittest -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi
echo "                 ...done"

echo "python3 -m unittest tests/XML/_unittest.py -v"
python3 -m unittest tests/XML/_unittest.py -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi
echo "                 ...done"

exit $ERR
