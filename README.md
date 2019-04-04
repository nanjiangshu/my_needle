# my_needle

## Description:
Global pairwise sequence alignment derived from the needle program in the EMBOSS package
Gapopens can be supplied under the sequence and enclosed in { }.
In addition, it allows pairwise alignment of a large number of sequences.
The software is developed in C/C++.

## Author

Nanjiang Shu

Senior researcher at SciLifeLab Stockholm

nanjiang.shu@scilifelab.se

## Usage:
```
my_needle [options] seqfile1 seqfile2

OPTIONS:
<pre>
      -list    FILE    Supply the pair list, one line a pair
      -mode     INT    Pair mode, (default: 0)
                       0: pairs are two filenames of sequences
                       1: pairs are two seqIDs of sequences
                       2: do all-to-all pairwise alignment given a fasta file
                          with multiple sequences
      -seqdb    STR    Indexed sequence database, this must be supplied when mode = 1
      -m        INT    Output alignment format, (default: 0)
                       0: full pairwise alignment in EMBOSS needle format
                       1: full pairwise alignment in Fasta format
      -table FILE Output tab delimited alignment info table to FILE
      -type     INT    Set the alignment type. (default: 1)
                       0: for dna, 
                       1: for aaseq
                       2: for shape strings
      -matrix  FILE    Supply substitution matrix file. default is
                       NUC.4.4:     for DNA alignment
                       BLOSUM62:    for amino acid alignment
                       ShagMatrix6: for shape string alignment
      -title1   STR    Force the title of the first seqfile
      -title2   STR    Force the title of the second seqfile
      -G      FLOAT    Penalty to open a gap, (default: 10.0)
      -E      FLOAT    Penalty to extend a gap, default: 0.5)
      -o       FILE    Output the alignment to outfile, (default: stdout)
      -eg     FLOAT    Set endgapopen, (default: 10.0)
      -ee     FLOAT    Set endgapextend, (default: 0.5)
      -debug   FILE    Print the debug information to FILE
      -endweight       Enable endweight, default= No
      -nogoarray       Do not use the gapopen array even supplied
      -threshold-seqdb-size INT
                       Threshold of seqdb file size in MB. (default: 1024)
                       If seqdb filesize < threshold, the whole file will be read
                       in to memory.
      -max-align-size INT
                       Maximum align size in MB. (default: 1024)
                       max_align_size = seqlength1 * seqlengh2
      -selpair INT     Pair selection method for mode 2, (default: 1)
                       0: sequentially 
                       1: randomly
      -wordsize INT    Wordsize for Kmer comparison, (default: 3)
      -pkmer|-printkmerscore yes|no
                       Whether print KMer score, (default: no)
     
      -h | --help      Print this help message and exit
</pre>

Note that the sequence file should be in FASTA format.

```

## Examples (in the "test" folder):

Align the first sequence in oneseq.fa to all sequences in the 10seq.fa

    $ ../my_needle oneseq.fa 10seq.fa

Align a large number of pairs of sequences by supplying sequence database and pairs of seqids

    $ ../my_needle -list PF00008.pairlist -mode 1 -seqdb PF00008.fa -m 1 -o outfile.txt


All to all pairwise alignment given a file with multiple sequences

    $ ../my_needle t2.fa -mode 2 -m 1 -table table1.txt


## Installation
First get the software by 

    $ git clone https://github.com/nanjiangshu/my_needle

Then change to the directory "my_needle" by

    $ cd my_needle

Then compile the software by running (note that you need gcc in order to
compile the software)

    $ make 

For making the debug version, type

    $ make debug


## Format of the input file

### Enhanced Fasta format

Gap penalties can be set under the sequence and enclosed in {}. A tag

    *gpo:  gap open
    *gpe:  gap extension
    *tgpe: terminal gap extension

can be set directly after the left bracket {. If no tag is set, the default is the gap open array.
A leading '#' will comment out the gap penalty array.

An example of input file with supplied gap open penalties:

```
>seq1
ASNLSKLFLSDSDA
{gpo: 10 -1000 10 10 10 10 10 10 10 10 10 10 10 10 }
{#gpe: 20 -1000 -10000 23 24 25 26 27 28 29 20 22 22 23 }
{#tgpe: 40 44 44 43 44 45 46 47 48 49 40 44 44 43 }
```


## Contact
Email: nanjiang.shu@scilifelab.se
