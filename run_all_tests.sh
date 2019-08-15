#!/bin/bash
ERR=0

echo "python3 -m unittest tests.POSCARloader._unittest -v"
python3 -m unittest tests.POSCARloader._unittest -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "python3 -m unittest tests.NearestNeighborsGenerator._unittest -v"
python3 -m unittest tests.NearestNeighborsGenerator._unittest -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "python3 -m unittest tests.argv._unittest -v"
python3 -m unittest tests.argv._unittest -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "python3 -m unittest tests/XML/_unittest.py -v"
python3 -m unittest tests/XML/_unittest.py -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "python3 -m unittest tests/JorGpi.py -v"
python3 -m unittest tests/JorGpi.py -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "python3 -m unittest tests/pickup.py -v"
python3 -m unittest tests/pickup.py -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

exit $ERR
