#!/bin/bash

wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz
tar -zxvf gsl-2.6.tar.gz
cd gsl-2.6
mkdir /home/yourname/gsl
./configure --prefix=/home/yourname/gsl
make -j2
make check
make install
