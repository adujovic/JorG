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
    ISINBASHRC=$(grep 'JORGPI_ASA_SRC' ${HOME}/.bashrc)
    if [ -z "${ISINBASHRC}" ]; then
      echo "export JORGPI_ASA_SRC=${PREFIX}/src" >> ${HOME}/.bashrc
    fi
else
    sudo mkdir -p /usr/src
    sudo cp -r JorGpi/asa /usr/src
    ISINBASHRC=$(sudo grep 'JORGPI_ASA_SRC' /etc/bash.bashrc)
    if [ -z "${ISINBASHRC}" ]; then
      echo "export JORGPI_ASA_SRC=/usr/src" | sudo tee -a /etc/bash.bashrc > /dev/null
    fi
fi
