#!/bin/bash

if [ -n "$JORGPI_ASA_SRC" ]; then
    echo "Simmulated annealing sources found in ${JORGPI_ASA_SRC}."
    exit 0
fi

while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
  --prefix )
    shift; PREFIX=$1
    ;;
  -v | --verbose )
    verbose=1
    ;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

if [ -n "$PREFIX" ]; then
    mkdir -p ${PREFIX}/src
    cp -r JorGpi/asa ${PREFIX}/src
    echo "export JORGPI_ASA_SRC=${PREFIX}/src" >> ${HOME}/.bashrc
else
    sudo mkdir -p /usr/src
    sudo cp -r JorGpi/asa /usr/src
    sudo echo "export JORGPI_ASA_SRC=/usr/src" >> /etc/bash.bashrc
fi
