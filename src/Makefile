CC = g++
CFLAGS = -std=c++11 -Wall -Wextra -g
TESTGEN = cxxtestgen
TESTFLAGS = --error-printer
TESTS = MyTest.h
TARGET = main
OBJS = main.o Utils.o Poly.o Basis.o Density.o

#Compiler le code
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -larmadillo

main.o: main.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp %.h
	$(CC) $(CFLAGS) -c $< -o $@

#Compiler les tests unitaires
test: tests
	./tests

tests.cpp: $(TESTS)
	$(TESTGEN) $(TESTFLAGS) -o $@ $(TESTS)

tests.o: tests.cpp
	$(CC) $(CFLAGS) -c $^ -o $@

tests: tests.o Poly.o Basis.o Utils.o Density.o
	$(CC) $(CFLAGS) $^ -o $@ -larmadillo

#Vérifier les fuites mémoires
memcheck: tests
	valgrind ./tests

#Préparation des données
pre-data:
	python rVals3d.py

#Tracer des graphes
plot2d: plot2d.txt
	python plot2d.py

plotpsi: Psi000.txt rVals_psi000.txt zVals_psi000.txt
	python Psi.py

.PHONY: clean

clean:
	rm -f *.o
	rm -f $(TARGET)
	rm -f tests.cpp
	rm -f tests
	rm -f plot.txt
	rm -f rVals.txt
	rm -f zVals.txt
	rm -f plot2d.txt
	rm -f plot3d.raw
	rm -f rVals3d.txt
	rm -f Psi*.txt
	rm -f rVals_psi*.txt
	rm -f zVals_psi*.txt



