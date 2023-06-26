all: prog 

prog: main.o hmmComparison.o
	g++ -std=c++11 main.o hmmComparison.o -o prog

main.o: main.cpp
	g++ -std=c++11 -c main.cpp

hmmComparison.o: hmmComparison.cpp
	g++ -std=c++11 -c hmmComparison.cpp

clean:
	rm *.o prog