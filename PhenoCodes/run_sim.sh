#!/bin/bash
#Este código ejecuta la simulación y crea las muestras para cada ParamCard. Debe encontrarse en la misma carpeta donde están las carpetas de las simulaciones. Ejecutar luego de name_sim.sh.
#Se realiza un ciclo para las ParamCards. Se debe escribir la dirección donde están ubicadas.
for i in $(ls /home/jgodoy/PhenoThesisProject/ParamCards); do
#Se chequea que se haya creado la carpeta de la ParamCard
if test -e ${i/.dat/}; then
#El código crea un archivo de comandos para madgraph
#Ejecuta la simulación
echo 'launch' ${i/.dat/} > mgcode_$COUNTER
#Usar Pythia y CMS
echo '1' >> mgcode_$COUNTER
echo '2' >> mgcode_$COUNTER
#Continuar
echo '0' >> mgcode_$COUNTER
echo '0' >> mgcode_$COUNTER
#Ejecutar MadGraph. se debe escribir la dirección donde el ejecutable de madgraph se encuentr ubicado
./home/jgodoy/PhenoThesisProject/Programs/MG5_aMC_v2_2_3/bin/mg5_aMC mgcode_$COUNTER &
else
#Si no se ha creado la simulación, se envía n mensaje.
echo Para el ParamCard $i no se ha creado la simulación
fi

done
