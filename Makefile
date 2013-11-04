.PHONY: clean main

VERSION=0.1-beta
CFLAGS=-O3

main:
	$(CXX) $(CFLAGS) -Wall -DVERSION=\"$(VERSION)\" -std=c++98 AmpliconClipper.cpp -o AmpliconClipper

clean:
	rm -f AmpliconClipper *.o
