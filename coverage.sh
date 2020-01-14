#!/bin/bash

(
cd JorGpi/tests
ln -s ../asa ./
ln -s ../_INPUT ./
ln -s ../_VASP ./
)

nosetests-3.4 --with-coverage --cover-xml -w JorGpi/tests
python-codacy-coverage -r JorGpi/tests/coverage.xml
rm JorGpi/tests/coverage.xml

cd JorGpi/tests
rm asa _INPUT _VASP
rm -r output
