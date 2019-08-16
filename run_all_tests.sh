#!/bin/bash
ERR=0

cd tests
ln -s ../asa
ln -s ../_INPUT
ln -s ../_VASP
echo "testing module POSCARloader"
python3 -m unittest test_POSCARloader -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "testing module NearestNeighborsGenerator"
python3 -m unittest test_generator -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "testing module argv"
python3 -m unittest test_argv -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "testing module xml"
python3 -m unittest test_xml -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "testing module JorGpi"
python3 -m unittest test_jorgpi -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "testing module pickup"
python3 -m unittest test_pickup -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

rm asa _INPUT _VASP
rm -r output
cd ..
exit $ERR
