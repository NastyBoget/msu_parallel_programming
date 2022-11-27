# module load SpectrumMPI
submit-polus-one:
	/usr/local/cuda/bin/nvcc -O3 -std=c++11 -arch=sm_35 --compiler-bindir /usr/bin/g++ -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm main.cu -o task4
	bsub -n 1 -gpu "num=2" -R "span[ptile=2]" -W 00:10 -o out.txt -e error.txt mpiexec ./task4 128 1 out\_1\_128\_1.txt

submit-polus:
	/usr/local/cuda/bin/nvcc -O3 -std=c++11 -arch=sm_35 --compiler-bindir /usr/bin/g++ -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm main.cu -o task4
	for N in 128 256 512 ; do \
		for p in 1 2 4 6 ; do \
			for i in {1..5} ; do \
				bsub -n $$p -gpu "num=2" -R "span[ptile=2]" -W 00:10 -o /dev/null -e /dev/null mpiexec ./task4 $$N 1 out\_$$p\_$$N\_1.txt ; \
				bsub -n $$p -gpu "num=2" -R "span[ptile=2]" -W 00:10 -o /dev/null -e /dev/null mpiexec ./task4 $$N out\_$$p\_$$N\_pi.txt ; \
			done \
		done \
	done
