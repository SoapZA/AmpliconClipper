.PHONY: clean main

main:
	g++ -O3 -std=c++0x AmpliconClipper.cpp -o AmpliconClipper

clean:
	rm -f AmpliconClipper *.o
