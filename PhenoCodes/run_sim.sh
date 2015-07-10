#!/bin/bash
#Este código crea las condiciones para la simulación en MadGraph. 
#Se crea un contador de prueba para verificar que se creen todas las simulaciones. Se ubica en la carpeta donde deben quedar las carpetas de las simulaciones.
COUNTER=0
echo El contador es $COUNTER
echo "import model mssm" > mgcode
echo "generate p p > ta1+ ta1- QCD=1" >> mgcode
echo "output auto" >> mgcode
echo "exit" >> mgcode

#Se realiza un ciclo que recorre todos las ParamCards. Acá se debe escribir la dirección donde se encuentran las ParamCards
for i in $(ls /home/jgodoy/PhenoThesisProject/ParamCards); do
#Este comando ejecuta MadGraph. Se debe escribir la dirección de el ejecutable de MadGraph. El tomando toma como entrada un archivo "mgcode" que debe ser modificado de acuerdo a la colisión que se quiera emular. El archivo importa el modelo "import mssm", luego genera la simulación (ej: generate p p > ta1+ ta1-), crea un output ("output auto") y luego sale de Madgraph ("exit"). El símbolo "&" ejecuta madgraph en paralelo.
./../../Programs/MG5_aMC_v2_2_3/bin/mg5_aMC /home/jgodoy/PhenoThesisProject/Simulations/mg_files/mgcode &
let COUNTER=COUNTER+1
#Se imprime en pantalla las simulaciones creadas
echo El contador es $COUNTER
done 
wait

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

#Este código ejecuta la simulación y crea las muestras para cada ParamCard. Debe encontrarse en la misma carpeta donde están las carpetas de las simulaciones. Ejecutar luego de name_sim.sh.
#Se realiza un ciclo para las ParamCards. Se debe escribir la dirección donde están ubicadas.
COUNTER=1
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
./../../Programs/MG5_aMC_v2_2_3/bin/mg5_aMC mgcode_$COUNTER &
else
#Si no se ha creado la simulación, se envía n mensaje.
echo Para el ParamCard $i no se ha creado la simulación
fi
let COUNTER=$COUNTER+1
done
wait
echo "Se ejecutaron las simulaciones"

#Se borra los archivos de ejecución. Ejecutar luego de run_sim.sh
rm mgcode*
rm py*
