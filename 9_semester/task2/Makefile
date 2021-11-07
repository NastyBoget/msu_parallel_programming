SRC=task2.cpp
CXX=mpic++
CXX_FLAGS=-std=c++11 -Wall -O2
OUT=sort


all: $(SRC) 
	$(CXX) $(SRC) -o $(OUT) $(CXX_FLAGS)

run: $(OUT)
	mpirun -n 4 ./$(OUT) 1000000

clean:
	rm -rf $(OUT)
