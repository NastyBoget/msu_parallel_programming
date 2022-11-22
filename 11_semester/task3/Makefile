submit-polus-one:
	g++ -O3 -std=c++11 -fopenmp main.cpp -o task3
	bsub -n 1 -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 ./task3 128 1 out\_1\_128\_1.txt
	g++ -O3 -std=c++11 -fopenmp main_mpi.cpp -o task3 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm
	bsub -n 4 -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./task3 128 1 out\_4\_128\_1.txt

# module load SpectrumMPI/10.1.0
submit-polus-consecutive:
	g++ -O3 -std=c++11 -fopenmp main.cpp -o task3
	for N in 128 256 512 ; do \
		for i in {1..5} ; do \
			bsub -n 1 -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 ./task3 $$N 1 out\_1\_$$N\_1.txt ; \
            bsub -n 1 -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 ./task3 $$N out\_1\_$$N\_pi.txt ; \
		done \
	done

submit-polus-parallel:
	g++ -O3 -std=c++11 -fopenmp main_mpi.cpp -o task3 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm
	for N in 128 256 512 ; do \
		for p in 1 4 8 16 32 ; do \
			for i in {1..5} ; do \
			  	bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./task3 $$N 1 out\_$$p\_$$N\_1.txt ; \
			  	bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./task3 $$N out\_$$p\_$$N\_pi.txt ; \
			done \
		done \
	done

submit-polus-consecutive-omp:
	g++ -O3 -std=c++11 -fopenmp main.cpp -o task3
	for N in 128 256 512 ; do \
		for i in {1..5} ; do \
			bsub -n 1 -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=4 ./task3 $$N 1 out\_1\_$$N\_1.txt ; \
            bsub -n 1 -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=4 ./task3 $$N out\_1\_$$N\_pi.txt ; \
		done \
	done

submit-polus-parallel-omp:
	g++ -O3 -std=c++11 -fopenmp main_mpi.cpp -o task3 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm
	for N in 128 256 512 ; do \
		for p in 1 4 8 16 32 ; do \
			for i in {1..5} ; do \
				bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=4 mpiexec ./task3 $$N 1 out\_$$p\_$$N\_1.txt ; \
                bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=4 mpiexec ./task3 $$N out\_$$p\_$$N\_pi.txt ; \
			done \
		done \
	done
