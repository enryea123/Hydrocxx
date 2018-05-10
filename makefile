INCS = `root-config --cflags`
LIBS = `root-config --libs`

comp: main.o vector.o density.o equation.o rungekutta.o shooting.o
	g++ -o program main.o vector.o density.o equation.o rungekutta.o shooting.o $(INCS) $(LIBS)

main.o: main.cxx vector.h
	g++ -o main.o -c main.cxx $(INCS)

vector.o: vector.cxx vector.h
	g++ -o vector.o -c vector.cxx

density.o: density.cxx density.h
	g++ -o density.o -c density.cxx

equation.o: equation.cxx equation.h
	g++ -o equation.o -c equation.cxx

rungekutta.o: rungekutta.cxx rungekutta.h
	g++ -o rungekutta.o -c rungekutta.cxx

shooting.o: shooting.cxx shooting.h
	g++ -o shooting.o -c shooting.cxx
