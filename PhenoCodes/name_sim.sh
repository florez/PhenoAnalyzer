#!/bin/bash
#Ese código toma los datos de las ParamCards y los integra en las muestras, para realizar la simulación. Además, la carpeta de cada simulación se le asigna el nombre de su respectivo ParamCard.
#Se crea un contador para verificar las ParamCards modificadas
COUNTER2=0
#Se realiza un ciclo en las ParamCards. Se debe escribir la ubicación de las ParamCards
for j in $(ls /home/jgodoy/PhenoThesisProject/ParamCards); do
#Se integran los datos de cada ParamCard para la simulación. Primero se escribe la dirección de cada ParamCard teniendo en cuenta que es un ciclo "$j". Luego se envía esos datos a la carpeta creada para esa simulación luego de haber ejecutado "create_sim.sh"
cat /home/jgodoy/PhenoThesisProject/ParamCards/$j > PROC_mssm_$COUNTER2/Cards/param_card.dat
#Se cambia el nombre de la carpeta para que concuerde con su respectivo ParamCard
mv PROC_mssm_$COUNTER2 ${j/.dat/}
#se imprime las muestras terminadas.
echo EL contador para cambiar es $COUNTER2
let COUNTER2=$COUNTER2+1
done

