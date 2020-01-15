#!/bin/bash

CMD=$(command -v gsl-config)
if [ -n "$CMD" ]; then
    echo "GSL found."
    exit 0
fi

for g in gsl-*; do
    if [ ! -d $g ]; then
        GSL_DIRECTORY=$g
        break
    fi
done

if [ -z "$GSL_DIRECOTRY" ]; then
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
    if [ ! -s gsl-latest.tar.gz ]; then
        echo "ERROR DOWNLOADING GSL!"
        exit 255
    fi
    tar -xzf gsl-latest.tar.gz
    rm -f gsl-latest.tar.gz
    $GSL_DIRECTORY=$(ls -d gsl*)
fi

while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
  --version )
    echo ${GSL_DIRECTORY:4:4}
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
  cd $GSL_DIRECTORY
  if [ -z $PREFIX ]; then
    ./configure "${COPTIONS}"
    make -j${NCPUS}
    sudo make install
    ISINBASHRC=$(sudo grep '/usr/local/lib' /etc/bash.bashrc)
    if [ -z "$ISINBASHRC" ]; then
      sudo echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib" >> /etc/bash.bashrc
    fi
  else
    ./configure --prefix=${PREFIX} "${COPTIONS}"
    make -j${NCPUS}
    make install
    ISINBASHRC=$(grep "${PREFIX}/lib" >> ${HOME}/.bashrc)
    if [ -z "$ISINBASHRC" ]; then
     echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:${PREFIX}/lib" >> ${HOME}/.bashrc
    fi
  fi  
  make clean
)
