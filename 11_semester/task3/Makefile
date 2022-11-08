all: main_mpi.cpp
	mpixlC -O3 -std=c++11 main_mpi.cpp -o task3

submit-polus-one:
	mpixlC -O3 -std=c++11 main_mpi.cpp -o task3
	mpisubmit.pl -p 10 -w 00:10 --stdout out.txt --stderr err.txt ./task3 -- 64 1 out\_10\_64\_1.txt ; \

# module load SpectrumMPI
submit-polus:
	mpixlC -O3 -std=c++11 main_mpi.cpp -o task3
	for N in 128 256 512 ; do \
			for p in 1 10 20 40 ; do \
				for i in {1..5} ; do \
				    mpisubmit.pl -p $$p -w 00:10 --stdout /dev/null --stderr /dev/null ./task3 -- $$N 1 out\_$$p\_$$N\_1.txt ; \
					mpisubmit.pl -p $$p -w 00:10 --stdout /dev/null --stderr /dev/null ./task3 -- $$N out\_$$p\_$$N\_pi.txt ; \
				done \
			done \
		done
