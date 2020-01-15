#!/bin/bash
ERR=0
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
  -c | --cpp | --full )
    shift; RUNCPP=1
    ;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

cd JorGpi/tests
ln -s ../asa ./
ln -s ../_INPUT ./
ln -s ../_VASP ./
echo "testing module POSCARloader"
python3 -m unittest test_poscar_loader -v
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

if [ -n "$RUNCPP" ]; then
  echo "testing module JorGpi"
  python3 -m unittest test_jorgpi -v
  TST=$?
  if [ "$TST" -ne "0" ]; then
      ERR=$TST
  fi
fi

echo "testing module pickup"
python3 -m unittest test_pickup -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

echo "testing module KPOINTS"
python3 -m unittest test_kpoints -v
TST=$?
if [ "$TST" -ne "0" ]; then
    ERR=$TST
fi

rm asa _INPUT _VASP
if [ -n "$RUNCPP" ]; then
  rm -r output
fi
exit $ERR
