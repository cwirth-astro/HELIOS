
NAME=test
MIXFILE=vertical_mix

# run the iteration for a sufficient number of iterations (e.g., 10)
for i in {3..10..1}
do
	# run HELIOS first
	python3 ./helios.py	-parameter_file=param_trappist.dat \
				-coupling_iteration_step $i

	# stops iteration after convergence is found
	if (( $i > 0 )) 
	then
		STOP=$(<helios_main_dir/output/${NAME}/${NAME}_coupling_convergence.dat)
   		echo -e "--> Converged? ${STOP} (1 = yes, 0 = no)"	
   		if ((${STOP}==1))
   		then
   			break
   		fi
	fi

	# run here your photochemical kinetics code
	/home/cwirth/FastChem/fastchem input/chemistry/fastchem_input/config.dat
done
