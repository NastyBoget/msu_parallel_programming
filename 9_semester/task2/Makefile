SRC=task2.cpp
CXX=mpicxx
CXX_FLAGS=-O3
OUT=sort


all: $(SRC) 
	$(CXX) $(SRC) -o $(OUT) $(CXX_FLAGS)

run: $(OUT)
	mpirun -n 4 ./$(OUT) 32 32

clean:
	rm -rf $(OUT)
