#!/bin/bash

#wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
#tar -xzf gsl-latest.tar.gz
#rm -f gsl-latest.tar.gz

CMD=$(command -v gsl-config)
if [ -n "$CMD" ]; then
    echo "GSL found."
    exit 0
fi
DIR=$(ls -d gsl*)
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
  --version )
    echo ${DIR:4:4}
    exit
    ;;
  --prefix )
    shift; PREFIX=$1
    ;;
  --coptions )
    shift; COPTIONS=$1
    ;;
  -v | --verbose )
    verbose=1
    ;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

NCPUS=$(grep -c proc /proc/cpuinfo)
if [ "$verbose" -eq "1" ]; then
    echo "I'm running $NCPUS procs"
fi
(
  cd $DIR
  ./configure --prefix=${PREFIX} "${COPTIONS}"
  make -j${NCPUS}
  if [ -z $PREFIX ]; then
    sudo make install
  else
    make install
  fi  
  make clean
)
