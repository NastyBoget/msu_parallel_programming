SRC=task3.cpp utils.cpp batcher_sort.cpp
CXX=mpixlcxx_r
OUT=task3


build: $(SRC)
	$(CXX) $(SRC) -o $(OUT) -O3 -std=c++11

run-single: build
	mpisubmit.bg -n 32 -w 02:00:00 --stdout result\_32\_$$K\_64x64\_net$$NET.txt --stderr result\_32\_$$K\_64x64\_net$$NET.err ./$(OUT) -- $$K 64 64 $$NET result\_32\_$$K\_64x64\_net$$NET\_out.txt

# make run-bluegene K=16,31,64 NET=1,2

run-bluegene: build
	for N in 128 256 512; do \
    	mpisubmit.bg -n 64 -w 00:30:00 --stdout result\_64\_$$K\_$$N\_net$$NET.txt --stderr result\_64\_$$K\_$$N\_net$$NET.err ./$(OUT) -- $$K $$N $$N $$NET result\_64\_$$K\_$$N\_net$$NET\_out.txt; \
		mpisubmit.bg -n 128 -w 00:10:00 --stdout result\_128\_$$K\_$$N\_net$$NET.txt --stderr result\_128\_$$K\_$$N\_net$$NET.err ./$(OUT) -- $$K $$N $$N $$NET result\_128\_$$K\_$$N\_net$$NET\_out.txt; \
		mpisubmit.bg -n 256 -w 00:10:00 --stdout result\_256\_$$K\_$$N\_net$$NET.txt --stderr result\_256\_$$K\_$$N\_net$$NET.err ./$(OUT) -- $$K $$N $$N $$NET result\_256\_$$K\_$$N\_net$$NET\_out.txt; \
		for P in 1 2 4 8 16 32; do \
        	mpisubmit.bg -n $$P -w 02:00:00 --stdout result\_$$P\_$$K\_$$N\_net$$NET.txt --stderr result\_$$P\_$$K\_$$N\_net$$NET.err ./$(OUT) -- $$K $$N $$N $$NET result\_$$P\_$$K\_$$N\_net$$NET\_out.txt; \
        done \
	done \

clean:
	rm -rf $(OUT) result\_* core*
