#!/bin/sh

rm $0

/prg/openmpi/bin/mpirun -np 8 "./main"

echo "

------------------
(program exited with code: $?)" 		


echo "Press return to continue"
#to be more compatible with shells like dash
dummy_var=""
read dummy_var
