#!/bin/bash
BF='\033[1m'
IT='\033[3m'
END='\033[0m'
DARKRED="\033[31m"
DARKGREEN="\033[32m"
DARKYELLOW="\033[33m"
DARKBLUE="\033[34m"

echo -e "$BF${DARKRED}Running tests:$END"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}-h :"
./argv_test.py -h
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}--incar INCAR --input POSCAR -N 3 -E Fe :"
./argv_test.py --incar INCAR --input POSCAR -N 3 -E Fe
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}--incar INCAR --input POSCAR -N 13 -E Ni :"
./argv_test.py --incar INCAR --input POSCAR -N 13 -E Ni
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}--incar INCAR --input POSCAR -N 3 -E Ni :"
./argv_test.py --incar INCAR --input POSCAR -N 3 -E Ni
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}--incar INCAR_tst1 --input POSCAR_tst1 -N 3 -E Fe :"
./argv_test.py --incar INCAR_tst1 --input POSCAR_tst1 -N 3 -E Fe
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}--symmetry -i POSCAR_Cs2F6Ni2  :"
./argv_test.py --symmetry -i POSCAR_Cs2F6Ni2 
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}-i POSCAR_tst1  :"
./argv_test.py -i POSCAR_tst1 
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}--incar INCAR_tst1 --input POSCAR_tst1 -N 2 -E Fe1 :"
./argv_test.py --incar INCAR_tst1 --input POSCAR_tst1 -N 2 -E Fe1
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}--incar _INCARs/INCAR_tst1 --input _POSCARs/POSCAR_tst1 -N 1 -E Fe1 :"
./argv_test.py --incar _INCARs/INCAR_tst1 --input _POSCARs/POSCAR_tst1 -N 1 -E Fe1
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}--incar _INCARs/INCAR_CsNiF --input _POSCARs/POSCAR_CsNiF -N 1 -E Ni -o output/J1 :"
./argv_test.py --incar _INCARs/INCAR_CsNiF --input _POSCARs/POSCAR_CsNiF -N 1 -E Ni -o output/J1
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}--symmetry -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF  :"
./argv_test.py --symmetry -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF 
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}-i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 6 -E Ni -o output/ASD --redundant :"
./argv_test.py -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 6 -E Ni -o output/ASD --redundant
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}-i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD --extra-dimentions \"2 2 2\" :"
./argv_test.py -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD --extra-dimentions "2 2 2"
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}-i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD -X \"2;2;2;2;4\" :"
./argv_test.py -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD -X "2;2;2;2;4"
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}-i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD -X \"2 ,2 ,2,4\" :"
./argv_test.py -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD -X "2 ,2 ,2,4"
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
echo -e "${IT}${DARKGREEN}  ./argv_test.py ${END}-i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD -X \"0 0 1\" :"
./argv_test.py -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD -X "0 0 1"
if [ "$?" == "0" ]; then
    echo -e "$BF$DARKBLUE Test succeed $END"
else
    echo -e "$BF$DARKYELLOW Test failed $END"
fi
echo "------------------------------------------------------------------------------------------------------------------------------------------------------"
