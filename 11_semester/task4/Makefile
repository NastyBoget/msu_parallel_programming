# module load OpenMPI
submit-polus-one:
	/usr/local/cuda/bin/nvcc -O3 -std=c++11 -arch=sm_35 --compiler-bindir /usr/bin/g++ -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm main.cu -o task4
	bsub -n 1 -gpu "num=2" -R "span[ptile=2]" -W 00:10 -o out.txt -e error.txt OMP_NUM_THREADS=1 mpiexec ./task4 128 1 out\_1\_128\_1.txt
