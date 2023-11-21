.PHONY: run

run :
	g++ -std=c++14 -Drestrict=__restrict__ -Wall -fdiagnostics-color=always -g -I./include ./src/*.cpp ./src/*.c -o launch
	./launch
