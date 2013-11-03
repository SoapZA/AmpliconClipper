.PHONY: clean main

CFLAGS=-O3

main:
	$(CXX) $(CFLAGS) -Wall -std=c++0x AmpliconClipper.cpp -o AmpliconClipper

clean:
	rm -f AmpliconClipper *.o
