#!/bin/bash

cd tests
ln -s ../asa ./
ln -s ../_INPUT ./
ln -s ../_VASP ./
cd ..

nosetests-3.4 --with-coverage --cover-xml -w tests
python-codacy-coverage -r tests/coverage.xml
rm tests/coverage.xml

cd tests
rm asa _INPUT _VASP
rm -r output
