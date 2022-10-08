all: task2.cpp
	mpixlC -O3 -std=c++11 task2.cpp -o task2

# module load SpectrumMPI
submit-polus:
	mpixlC -O3 -std=c++11 task2.cpp -o task2

	for eps in 3e-5 5e-6 1.5e-6 ; do \
		for p in 1 4 16 64 ; do \
			for i in {1..10} ; do \
				mpisubmit.pl -p $$p -w 00:10 --stdout /dev/null --stderr /dev/null ./task2 -- $$eps out\_$$p\_$$eps.txt ; \
			done \
		done \
	done

submit-bluegene:
	mpixlcxx -O3 task2.cpp -o task2

	for eps in 1e-4 2e-5 8e-6 ; do \
		for p in 1 4 16 64 ; do \
		  	for i in {1..10} ; do \
				mpisubmit.bg -n $$p -w 00:10:00 --stdout /dev/null --stderr /dev/null ./task2 -- $$eps out\_$$p\_$$eps.txt ; \
			done \
		done \
	done
