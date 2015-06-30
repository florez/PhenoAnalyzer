#!/bin/bash
#Este archivo crea las nuplas de Delphes a partir de las muestras de MadGraph creadas. Este código debe estar ubicado en la carpeta donde se quiere tener las nutplas.

prueba1=0
#Contador de las ParamCards
cont1=0
#Ciclo en las ParamCards, escribir la carpeta donde se encuentran
for i in $(ls /home/jgodoy/PhenoThesisProject/ParamCards); do
        #Variable con el nombre de la ParamCard
	archivo=${i/.dat/}
        #Verifica que exista la simulación para esa ParamCard. Escribir la carpeta donde están las carpetas de las muestras realizadas.
	if test -e /home/jgodoy/PhenoThesisProject/Simulations/mg_files/$archivo; then
                #Crear carpeta con el nombre de la ParamCard para guardar las ntuplas correspondientes a la misma
		mkdir $archivo
                #Contador de las corridas de cada ParamCard (*run_01...)
		cont2=1
                #Prueba de la existencia de más corridas. 0 Verdadero, 1 Falso.
		prueba2=0
                #Ciclo para la s corridas de dicho ParamCard. Mientras existan más corridas
		while [ "$prueba2" -eq "0" ]; do 
			if [ "$cont2" -lt "10" ]; then
                                #Comprueba que exista el archivo .hep para el ParamCard y el run en el ciclo. La ruta debe seguir la carpeta donde están las carpetas de las muestras de MadGraph creadas.
				if test -e /home/jgodoy/PhenoThesisProject/Simulations/mg_files/$archivo/Events/run_0$cont2/tag_1_pythia_events.hep*; then
					#Prueba exitosa en la consola
                                        echo existe run $cont2 para $archivo
					#Descomprime el archivo .hep
					gunzip /home/jgodoy/PhenoThesisProject/Simulations/mg_files/$archivo/Events/run_0$cont2/tag_1_pythia_events.hep.gz
					#Copia el archivo .hep descomprimido a la carpeta creada anteriormente donde se ubicará la ntupla para dicha ParamCard. Primero va la dirección con el archivo .hep y luego la dirección de la carpeta creada.
					cp /mg_files/$archivo/Events/run_0$cont2/tag_1_pythia_events.hep /home/jgodoy/PhenoThesisProject/Simulations/ntuple_delphes/$archivo
					temp="$archivo"_run_"$cont2"
					#Ejecuta Delphes y crea la ntupla dentro de la carpeta de la respectiva ParamCard. Se le asigna como nombre a la tupla el mismo nombre del archivo .hep con la terminación .root 
					./home/jgodoy/PhenoThesisProject/Programs/Delphes-3.2.0/DelphesSTDHEP /home/jgodoy/PhenoThesisProject/Programs/Delphes-3.2.0/cards/delphes_card_CMS.tcl $temp.root $archivo/tag_1_pythia_events.hep &
					let cont2=$cont2+1
				else
					let prueba2=1
					echo cambia prueba 2 a falso
				fi
			else
				if test -e /home/jgodoy/PhenoThesisProject/Simulations/mg_files/$archivo/Events/run_$cont2/tag_1_pythia_events.hep*; then
                                        echo existe run $cont2 para $cont1
					gunzip /home/jgodoy/PhenoThesisProject/Simulations/mg_files/$archivo/Events/run_$cont2/tag_1_pythia_events.hep.gz
					cp /home/jgodoy/PhenoThesisProject/Simulations/mg_files/$archivo/Events/run_$cont2/tag_1_pythia_events.hep /home/jgodoy/PhenoThesisProject/Simulations/ntuple_delphes/$archivo
					temp="$archivo"+"_run_"+"$cont2"
					./home/jgodoy/PhenoThesisProject/Programs/Delphes-3.2.0/DelphesSTDHEP /home/jgodoy/PhenoThesisProject/Programs/Delphes-3.2.0/cards/delphes_card_CMS.tcl $temp.root $archivo/tag_1_pythia_events.hep &
                                
					let cont2=$cont2+1
                                else
                                        let prueba2=1
					echo prueba 2 falsa
				fi
			fi
		done 
		let cont1=$cont1+1
	else
		let prueba1=1
		echo prueba1 falsa
	fi
done

