# AmpliconClipper
AmpliconClipper performes post-alignment PCR Primer removal. In a multitude of situations, pre-alignment clipping is not feasible; either because the primers 
themselves become progressively shorter or the sequencing platform is rather noisy (AmpliconClipper was designed with Pacific Biosciences data in mind.)

## Features
 - Takes SAM as input
 - Can specify either primers or amplicons directly
 - Can specify whether forward resp. reverse reads can only come from specific amplicons
 - Can discard reads with a below fraction coverage of an amplicon

#### PREREQUISITES TO BUILD:
 - You will require a C++ compiler that can hanlde C++0x (which includes C++11 compliant compilers, obviously)

## Building:
```
git clone git://github.com/SoapZA/AmpliconClipper.git;
cd AmpliconClipper;
make;
```

## Running:
See the helptext:
```
AmpliconClipper
         by David Seifert 2013

Options:
         -i   : input SAM file
         -o   : output SAM file
         -a   : input amplicon file
         -S   : write statistics of reads
         -I   : maximum insertions threshold
         -D   : maximum deletions threshold
         -C   : minimum coverage threshold as a fraction of the amplicon insert
        --mD  : maximum adjacent deletions
        --ref : provide the reference genome in order to locate amplicons via their primers
```

## LICENSE:
GNU GPLv3 http://www.gnu.org/licenses/gpl-3.0
