# AmpliconClipper
AmpliconClipper performes post-alignment PCR Primer removal. In a multitude of situations, pre-alignment clipping is not feasible; either because the primers 
themselves become progressively shorter or the sequencing platform is rather noisy (AmpliconClipper was designed with Pacific Biosciences data in mind).

## Features
 - Takes SAM as input
 - Can specify either primers or amplicons directly
 - Can specify whether forward resp. reverse reads can only come from specific amplicons
 - Can discard reads with a below fraction coverage of an amplicon

#### PREREQUISITES TO BUILD:
 - You will require a C++ compiler that can handle C++98

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
         -S   : (optional) write statistics of reads
         -I   : (optional) maximum insertions threshold
         -D   : (optional) maximum deletions threshold
         -C   : (optional) minimum coverage threshold as a fraction of the amplicon insert
        --mD  : (optional) maximum adjacent deletions
        --ref : required if you provide a list of primers. Amplicons are found by locating primers on the reference genome.
```

#### Specifying amplicons:
In order to specify amplicons, you have two options:
 - Either provide a file with the start, end, start of the insert and end of the insert of the amplicon. Additionally, specify in which direction reads can 
originate from this amplicon, i.e. '+' for this amplicon only yielding reads in the forward (sense) direction and '-' for reads in the anti-sense direction. 
Provide '0' if both sense and anti-sense reads can originate from this amplicon. The last column can contain a capital 'X', such that reads from this amplicon 
are discarded.
 - You can also give AmpliconClipper a FASTA file with primers. Every FASTA ID should occur not more or less than twice and AmpliconClipper will locate all 
amplicon parameters by searching the reference genome.

## LICENSE:
GNU GPLv3 http://www.gnu.org/licenses/gpl-3.0
