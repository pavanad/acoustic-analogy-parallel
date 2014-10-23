
CC=g++
CCMACOS=g++-mp-4.7
FLAGS=-pedantic -Wall -lpthread -lm -fopenmp -DUse_Complex

all: compile

compile: src/*.cpp src/*.h
	$(CCMACOS) $(FLAGS) -o surf src/*.cpp

clean:
	rm surf
	rm output/*.*

run:
	./surf -t 8 -o 360 -serial -parallel -output input01
	./surf -t 8 -o 720 -serial -parallel -output input02
	./surf -t 8 -o 1440 -serial -parallel -output input03

.PHONY: all clean