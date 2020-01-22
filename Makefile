all: solve 

solve: symplex.o main.o
		g++ symplex.o main.o -o symplex -o main 

symplex.o: symplex.c 
		g++ -c symplex.c

main.o: main.c 
		g++ -c main.c



clean:
	del -rf *.o *.exe
