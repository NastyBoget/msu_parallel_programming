SRC=task2.cpp
CXX=mpixlcxx_r
OUT=task2


build: $(SRC)
	$(CXX) $(SRC) -o $(OUT) -O3 -std=c++11

run-single: build
	mpisubmit.bg -n 128 -w 00:10:00 --stdout result\_128\_4096.txt --stderr result\_128\_4096.err ./$(OUT) -- 4096 4096

run-bluegene: build
	for N in 1024 2048 4096 8192; do \
    	mpisubmit.bg -n 64 -w 00:30:00 --stdout result\_64\_$$N.txt --stderr result\_64\_$$N.err ./$(OUT) -- $$N $$N; \
		mpisubmit.bg -n 128 -w 00:10:00 --stdout result\_128\_$$N.txt --stderr result\_128\_$$N.err ./$(OUT) -- $$N $$N; \
		mpisubmit.bg -n 256 -w 00:10:00 --stdout result\_256\_$$N.txt --stderr result\_256\_$$N.err ./$(OUT) -- $$N $$N; \
		for P in 1 2 4 8 16 32; do \
        	mpisubmit.bg -n $$P -w 01:30:00 --stdout result\_$$P\_$$N.txt --stderr result\_$$P\_$$N.err ./$(OUT) -- $$N $$N ; \
        done \
	done \

unite-results:
	cat result\_*.txt > results.txt

clean:
	rm -rf $(OUT) result\_* core*
