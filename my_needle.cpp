/* @funcstatic embAlignGetScoreNWMatrix ***************************************
**
** Returns score of the optimal global or overlap alignment for
** the specified path matrix for Needleman Wunsch
**
** @param [r] ix [const float*] Gap scores array, ix(i,j) is the best score
**                              given that a(i) is aligned to a gap
**                              (in an insertion with respect to b)
** @param [r] iy [const float*] Gap scores array, iy(i,j) is the best score
**                              given that b(i) is in an insertion
**                              with respect to a
**
** @param [r] m [const float*] Match scores array, m(i,j) is the best score
**                             up to (i,j) given that a(i) is aligned to b(j)
** @param [r] lena [ajint] length of the first sequence
** @param [r] lenb [ajint] length of the second sequence
** @param [w] start1 [ajint *] start of alignment in first sequence
** @param [w] start2 [ajint *] start of alignment in second sequence
** @param [r] endweight [AjBool] whether the matrix was built for
**                               a global or overlap alignment
**
** @return [float] optimal score
******************************************************************************/
/*
 * =====================================================================================
 *       Filename:  my_needle.cpp
 *    Description:  global sequence alignment with position specific gap
 *                  penalties
 *        Created:  09/09/2010 11:26:51 PM CEST
 *       Compiler:  g++
 *         Author:  Nanjiang Shu, nanjiang.shu@scilifelab.se
 *        Company:  Department of Biochemistry and Biophysics, Stockholm Univesity
 * =====================================================================================
 */
#include <cstdio>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <cmath>

#include <string>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include <cassert>
#include <ctime>
#include <climits>
#include "array.h"
using namespace std;

typedef unsigned char BYTE;
struct Pair{/*{{{*/
    string mem1;
    string mem2;
};/*}}}*/
struct dbindex{/*{{{*/
    int dbfileindex;
    long offset;
    unsigned long size;
};/*}}}*/

/*the gapopenarray, 
 * if want to open a gap at position 4, set the gapopen score at 4th position
 * to be a low value
 * e.g. 
 *
 * ACDEGD
 * AVDEGD
 *
 * if you want an alignment 
 * ACD-EGD
 * AVDEGD-
 *
 * set the gapopen value of seq1 gapopen[3] = -1000
 *
 * */

#ifndef WHITE_SPACE
#define WHITE_SPACE " \t\r\n"
#endif

#ifndef SEQALIGN
#define SEQALIGN
/*for the alignment result*/
#define IDT 3 /* identity*/
#define SIM 2 /* similar, positive*/
#define MIS 1 /* mismatch*/
#define GAP 0 /* gap*/
#endif

#ifndef ERROR_CODE
#define READ_FILE_ERROR  1
#define WRITE_FILE_ERROR 2
//#define READ_OVERFLOW    5
//#define WRITE_OVERFLOW   5
//#define UNDERFLOW        10
#endif /*ERROR_CODE*/

#ifndef SEQ_TYPE
#define SEQ_TYPE
#define DNA_SEQ 0 // dna sequence
#define AA_SEQ  1 // amino acid sequence
#define SHAPE_SEQ 2 // shape string sequence
#define UNKNOWN_SEQ_TYPE -1
#endif

typedef int ajint;
typedef unsigned int ajuint;
typedef ajint AjBool;
typedef signed char int8;

#define ajFalse 0
#define ajTrue 1

#define U_FEPS 1.192e-6F         /* 1.0F + E_FEPS != 1.0F */
#define U_DEPS 2.22e-15          /* 1.0 +  E_DEPS != 1.0  */

#undef DIAG
#undef LEFT
#undef DOWN
#define DIAG 0
#define LEFT 1
#define DOWN 2

#undef CHAR_GAP
#define CHAR_GAP '-'

#define E_FPEQ(a,b,e) (((b - e) < a) && (a < (b + e)))

#undef SIZE_TITLE
#define SIZE_TITLE 100

#undef LONGEST_SEQ 
int LONGEST_SEQ=80000;
float gapopen = 10.0;
float gapextend = 0.5;
float endgapopen = 10.0;
float endgapextend = 10.0;
int   seq_type = AA_SEQ; /* sequence type, default is amino acid sequence*/
int   nchar = 50; /*maximum length of sequence to print at each line*/
bool isCheckSeqIDTClass = true;
bool isShowProgress = true;
bool isRandSelection = true;
int  method_select_pair = 1;

#define HIGH_COVERAGE 0
#define HIGH_PRECISION 1

int  KMerThresholdMode = HIGH_PRECISION;
bool isPrintKMerBitScore = false;

int wordsize = 3;

int maxSeqIDTClass[]={
    10000,  /*0  - 10*/
    10000,  /*10 - 20*/
    10000,  /*20 - 30*/
    10000,  /*30 - 40*/
    10000,  /*40 - 50*/
    10000,  /*50 - 60*/
    10000,  /*60 - 70*/
    10000,  /*70 - 80*/
    10000,  /*80 - 90*/
    10000   /*90 - 100*/
};
float binSeqIDTClass[] = {
     0.0, 10.0,
    10.0, 20.0,
    20.0, 30.0,
    30.0, 40.0,
    40.0, 50.0,
    50.0, 60.0,
    60.0, 70.0,
    70.0, 80.0,
    80.0, 90.0,
    90.0, 100.0
};
int numSeqIDTClass = sizeof(maxSeqIDTClass)/sizeof(int);

const int para20[]={
    1,
    20,
    20*20,
    20*20*20,
    20*20*20*20,
    20*20*20*20*20,
    20*20*20*20*20*20,
    20*20*20*20*20*20*20,
};

const int aacode[]=
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*first 65 chars*/
/*unique AA, ARNDCQEGHILKMFPSTWYVBZX */
/*A,B, C,D,E,F ,G,H,I,J, K, L, M, N,O, P, Q,R,S ,T, U, V, W, X, Y, Z */
  0,20,4,3,6,13,7,8,9,22,11,10,12,2,22,14,5,1,15,16,22,19,17,22,18,21
 -1,-1,-1,-1,-1,-1,-1, 
/*a,b, c,d,e,f ,g,h,i,j, k, l, m, n,o, p, q,r,s ,t, u, v, w, x, y, z */
  0,20,4,3,6,13,7,8,9,22,11,10,12,2,22,14,5,1,15,16,22,19,17,22,18,21
};
#define NUM_BLOSUM 24   //size of BLOSUM matrix
#define NUM_NUC 15
#define NUM_SHAG 12

int MAX_ALIGN_SIZE=1024;  /*1024 MB*/
int THRESHOLD_DB_FILESIZE = 0; /*0 MB*/

/* Experiments show that the most time consuming part is the alignment, not 
 * file reading. The gain by reading all sequences in to memory is very little.
 * so just use fseek and fread technique.
 *
$  /data3/program/my_needle/my_needle -mode 1 -l t3.pair -seqdb pfamfullseq_uniq.fasta -o t3.txt -threshold-seqdb-size 0 -m 1   
# Read Database Index with 1393148 records costs 4.51 seconds.
# Align 3000 pairs of sequences costs 9.04 seconds.

/data3/program/my_needle/my_needle -mode 1 -l t3.pair -seqdb pfamfullseq_uniq.fasta -o t3.1.txt -threshold-seqdb-size 1024 -m 1   

# GetIDSeqMap of 1393148 sequences costs 7.49 seconds.
# Align 3000 pairs of sequences costs 8.95 seconds.   
*/



struct AlignFactor/*{{{*/
{
	float score;
	float zScore;
	float pozScore;
	double eValue;
	int   idt_cnt;
	float identity;
	float identity_short;
	int   sim_cnt;
	float similarity;
	float similarity_short;
	int   gap_cnt;
	float gapPercent;
};/*}}}*/

const char ALPHABET_BLOSUM[]="ARNDCQEGHILKMFPSTWYVBZX*";
const char ALPHABET_NUC[]="ATGCSWRYKMBVHDN";
const char ALPHABET_SHAG[]="AKRSUVTGHBC-";

int outformat = 0; /*set by the -m option 
                0: full pairwise alignment 
                1: full pairwise alignment in Fasta format, they are placed sequentially */
int pairlistmode=0;/*set by the -mode option
                0: pairs are two filenames of sequences
                1: pairs are two seqIDs of sequences
                2: all-to-all alignment given a file  with multiple sequences
                */

const int blosum62[][NUM_BLOSUM] =/*{{{*/
{  // substitution matrix for amino acid sequences 
   /*  A    R  N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
    {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
    { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
    { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
    { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
    {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
    { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
    { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
    {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
    { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
    { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
    { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
    { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
    { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
    { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
    { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
    {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
    {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
    { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
    { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
    {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
    { -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
    { -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
    {  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
    { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
};
/*}}}*/
const int nuc44[][NUM_NUC] =/*{{{*/
{   // substitution matrix for DNA 
    //# This matrix was created by Todd Lowe   12/10/92
    //# Uses ambiguous nucleotide codes, probabilities rounded to nearest integer
    //# Lowest score = -4, Highest score = 5
    //A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
    { 5, -4 , -4 , -4 , -4 ,  1 ,  1 , -4 , -4 ,  1 , -4 , -1 , -1 , -1 , -2 },
    {-4,  5 , -4 , -4 , -4 ,  1 , -4 ,  1 ,  1 , -4 , -1 , -4 , -1 , -1 , -2 },
    {-4, -4 ,  5 , -4 ,  1 , -4 ,  1 , -4 ,  1 , -4 , -1 , -1 , -4 , -1 , -2 },
    {-4, -4 , -4 ,  5 ,  1 , -4 , -4 ,  1 , -4 ,  1 , -1 , -1 , -1 , -4 , -2 },
    {-4, -4 ,  1 ,  1 , -1 , -4 , -2 , -2 , -2 , -2 , -1 , -1 , -3 , -3 , -1 },
    { 1,  1 , -4 , -4 , -4 , -1 , -2 , -2 , -2 , -2 , -3 , -3 , -1 , -1 , -1 },
    { 1, -4 ,  1 , -4 , -2 , -2 , -1 , -4 , -2 , -2 , -3 , -1 , -3 , -1 , -1 },
    {-4,  1 , -4 ,  1 , -2 , -2 , -4 , -1 , -2 , -2 , -1 , -3 , -1 , -3 , -1 },
    {-4,  1 ,  1 , -4 , -2 , -2 , -2 , -2 , -1 , -4 , -1 , -3 , -3 , -1 , -1 },
    { 1, -4 , -4 ,  1 , -2 , -2 , -2 , -2 , -4 , -1 , -3 , -1 , -1 , -3 , -1 },
    {-4, -1 , -1 , -1 , -1 , -3 , -3 , -1 , -1 , -3 , -1 , -2 , -2 , -2 , -1 },
    {-1, -4 , -1 , -1 , -1 , -3 , -1 , -3 , -3 , -1 , -2 , -1 , -2 , -2 , -1 },
    {-1, -1 , -4 , -1 , -3 , -1 , -3 , -1 , -3 , -1 , -2 , -2 , -1 , -2 , -1 },
    {-1, -1 , -1 , -4 , -3 , -1 , -1 , -3 , -1 , -3 , -2 , -2 , -2 , -1 , -1 },
    {-2, -2 , -2 , -2 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 }
};/*}}}*/
const int shag6[][NUM_SHAG] =/*{{{*/
{   // substitution matrix for shape strings
//# substitution matrix in integer format, scale = 9
//# modified on 2007-03-18, 
//# the score of A-K should be positive and U-V as well
//# in that case, manually modify 
//# A-A = 3
//# A-K = 2
//# A-U/V = 1
//# K-U/V = 1
//# R-S   = 1
//# R-U/V = 1
//    A  K  R   S   U   V   T    G    H    B   C    - 
    { 3, 2,-4, -6,  1,  1, -2,  -3,   2,  -5, -3,  -1, },
    { 2, 5,-1, -3,  1,  1,  0,   0,   1,  -2,  1,   0, },
    {-4,-1, 4,  1,  1,  1, -1,   0,  -3,   2,  0,   0, },
    {-6,-3, 1,  4,  1,  1, -3,   0,  -6,   3, -1,  -1, },
    { 1, 1, 1,  1,  9,  5,  1,   1,  -2,   0,  4,   0, },
    { 1, 1, 1,  1,  5,  8,  1,   0,  -2,   1,  4,   0, },
    {-2, 0,-1, -3,  1,  1,  9,   2,  -2,  -2,  7,   0, },
    {-3, 0, 0,  0,  1,  0,  2,  11,  -2,   0,  6,   0, },
    { 2, 1,-3, -6, -2, -2, -2,  -2,   2,  -5, -2,  -1, },
    {-5,-2, 2,  3,  0,  1, -2,   0,  -5,   3, -1,   0, },
    {-3, 1, 0, -1,  4,  4,  7,   6,  -2,  -1,  6,   0, },
    {-1, 0, 0, -1,  0,  0,  0,   0,  -1,   0,  0,   0, }
};/*}}}*/

char alphabet[500]          = "";
AjBool endweight            = false;
bool show                   = false;
string  matrixFile = "";
int  sizeAlphabet           = 0;
bool isPrintTraceMatrix     = false;
bool noGapOpenArray         = false;
int alnType                 = AA_SEQ;

char usage[]/*{{{*/="\n\
usage: my_needle [options] seqfile1 seqfile2\n\
\n\
Pairwise sequence alignment derived from the needle program in the EMBOSS package\n\
Gapopens can be supplied under the sequence and enclosed in { }\n\
\n\
Options:\n\
 -list    FILE    Supply the pair list, one line a pair\n\
 -mode     INT    Pair mode, (default: 0)\n\
                  0: pairs are two filenames of sequences\n\
                  1: pairs are two seqIDs of sequences\n\
                  2: do all-to-all pairwise alignment given a fasta file\n\
                     with multiple sequences\n\
 -seqdb    STR    Indexed sequence database, this must be supplied when mode = 1\n\
 -m        INT    Output alignment format, (default: 0)\n\
                  0: full pairwise alignment in EMBOSS needle format\n\
                  1: full pairwise alignment in Fasta format\n\
 -table FILE Output tab delimited alignment info table to FILE\n\
 -type     INT    Set the alignment type. (default: 1)\n\
                  0: for dna, \n\
                  1: for aaseq\n\
                  2: for shape strings\n\
 -matrix  FILE    Supply substitution matrix file. default is\n\
                  NUC.4.4:     for DNA alignment\n\
                  BLOSUM62:    for amino acid alignment\n\
                  ShagMatrix6: for shape string alignment\n\
 -title1   STR    Force the title of the first seqfile\n\
 -title2   STR    Force the title of the second seqfile\n\
 -G      FLOAT    Penalty to open a gap, (default: 10.0)\n\
 -E      FLOAT    Penalty to extend a gap, default: 0.5)\n\
 -o       FILE    Output the alignment to outfile, (default: stdout)\n\
 -eg     FLOAT    Set endgapopen, (default: 10.0)\n\
 -ee     FLOAT    Set endgapextend, (default: 0.5)\n\
 -debug   FILE    Print the debug information to FILE\n\
 -endweight       Enable endweight, default= No\n\
 -nogoarray       Do not use the gapopen array even supplied\n\
 -threshold-seqdb-size INT\n\
                  Threshold of seqdb file size in MB. (default: 1024)\n\
                  If seqdb filesize < threshold, the whole file will be read\n\
                  in to memory.\n\
 -max-align-size INT\n\
                  Maximum align size in MB. (default: 1024)\n\
                  max_align_size = seqlength1 * seqlengh2\n\
 -selpair INT     Pair selection method for mode 2, (default: 1)\n\
                  0: sequentially \n\
                  1: randomly\n\
 -wordsize INT    Wordsize for Kmer comparison, (default: 3)\n\
 -pkmer|-printkmerscore yes|no\n\
                  Whether print KMer score, (default: no)\n\
\n\
 -h | --help      Print this help message and exit\n\
\n\
Created on 2010-09-10, updated 2012-05-30, Nanjiang\n\
\n\
Examples:\n\
    # align the first sequence in test1.fa to all sequences in the test2.fa\n\
    my_needle test1.fa test2.fa\n\
\n\
    # align a large number of pairs of sequences by supplying sequence db and\n\
    # pairs of seqids\n\
    my_needle -list pairlist.txt -mode 1 -seqdb seqdb -m 1 -o outfile.txt\n\
\n\
    # all to all pairwise alignment given a file with multiple sequences\n\
    my_needle t2.fa -mode 2 -m 1 -table table1.txt\n\
";
/*}}}*/

char explanation_table_format[]="\
# Seq1       Sequence 1 in alignment\n\
# Seq2       Sequence 2 in alignment\n\
# IDT0       numIDTRes / alnLength\n\
# SIM0       numSIMRes / alnLength\n\
# AlnLength  length of the alignment\n\
# Len1       length of sequence 1 in the pair\n\
# Len2       length of sequence 2 in the pair\n\
# Score      alignment score\n\
# N_IDT      number of identical residues\n\
# N_SIM      number of similar (including identical) residues\n\
# N_GAP      number of gaps\n\
# IDT1       numIDTRes / min(len1, len2)\n\
# IDT2       numIDTRes / (alnLength - N_GAP)\n\
#\n";
void PrintHelp() {/*{{{*/
    fprintf(stdout,"%s\n",usage);
}     /*}}}*/
void PrintNeedleAlignmentHeader(string &programname,struct tm* timeinfo,string& outfile, int argc, char **argv,FILE *fpout)/*{{{*/
{
    fprintf(fpout,"########################################\n");
    fprintf(fpout,"# Program: %s\n", programname.c_str());
    fprintf(fpout,"# Rundate: %s", asctime(timeinfo));
    fprintf(fpout,"# Commandline:  ");
    for (int jj = 0; jj < argc; jj++) {
        fprintf(fpout,"%s ", argv[jj]);
    }
    fprintf(fpout,"\n");
    fprintf(fpout,"# Align_format: %s\n", "srspair");
    fprintf(fpout,"# Report_file: %s\n", outfile.c_str());
    fprintf(fpout,"########################################\n");
    fprintf(fpout,"\n");
}/*}}}*/
bool IsFileExist(const char *name)/*{{{*/
{
    struct stat st;
    return (stat(name,&st) == 0); 
}/*}}}*/ 
bool IsInCharSet(const char ch, const char *charSet, int n /*= 0 */)/*{{{*/
/*****************************************************************************
 * check if the character "ch" is in charSet,
 ****************************************************************************/
{
	if(n == 0) {
        n = strlen(charSet);
    }
    int i;
	for(i = 0 ;i < n ; i ++) {
		if(ch == charSet[i]) {
			return true;
        }
	}
	return false;
}/*}}}*/
int GetFileSize(const char *name)/*{{{*/
{
    struct stat st;
    if (stat(name,&st) == 0){
        return st.st_size; 
    }else{
        cerr << "File " 
            << name
            << " does not exist." << endl;
        return -1;
    }
}/*}}}*/ 
string int2string(int number)/*{{{*/
{
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}/*}}}*/
void InitAlignFactor(AlignFactor *pAlignFactor)/*{{{*/
{
    pAlignFactor -> score            = 0.0;
    pAlignFactor -> zScore           = 0.0;
    pAlignFactor -> pozScore         = 0.0;
    pAlignFactor -> eValue           = 0.0;
    pAlignFactor -> idt_cnt          = 0;
    pAlignFactor -> identity         = 0.0;
    pAlignFactor -> identity_short   = 0.0;
    pAlignFactor -> sim_cnt          = 0;
    pAlignFactor -> similarity       = 0.0;
    pAlignFactor -> similarity_short = 0.0;
    pAlignFactor -> gap_cnt          = 0;
    pAlignFactor -> gapPercent       = 0.0;

}
/*}}}*/
bool IsAllGroupFilled(vector <int> &cntSeqIDTClass)/*{{{*/
{
    for (int i = 0 ; i < numSeqIDTClass ; i ++) {
        if (cntSeqIDTClass[i] < maxSeqIDTClass[i]){
            return false;
        }
    }
    return true;
}
/*}}}*/

int GetSeqIDTClass(float seqIdentity, float* binSeqIDTClass, int numSeqIDTClass){/*{{{*/
    for (int i = 0; i < numSeqIDTClass; i++){
        //printf("%g   %.3g  %g\n", binSeqIDTClass[i], seqIdentity, binSeqIDTClass[i+1]);
        if (seqIdentity >= binSeqIDTClass[2*i] && seqIdentity < binSeqIDTClass[2*i+1]){
            return i;
        }
    }
    return numSeqIDTClass;
}/*}}}*/
float GetThresholdKMerScore(float minimalSeqIdentity){/*{{{*/
    /*sequence pair with kMerBitScore < threshold_kmerScore are not likely to
     * have sequence identity > minimalSeqIdentity*/
    /*To be developed*/
    if (KMerThresholdMode == HIGH_COVERAGE ){
        if (minimalSeqIdentity >=95.0){
            return -0.3;
        } else if (minimalSeqIdentity>=90.0){
            return -0.5;
        } else if (minimalSeqIdentity>=80.0){
            return -1.0;
        }else if (minimalSeqIdentity>=70.0){
            return -1.3;
        }else if (minimalSeqIdentity>=60.0){
            return -1.5;
        }else if (minimalSeqIdentity>=50.0){
            return -1.6;
        }else if (minimalSeqIdentity>=40.0){
            return -1.8;
        }else if (minimalSeqIdentity>=30.0){
            return -2.0;
        } else if (minimalSeqIdentity>=20.0){
            return -2.2;
        }else if (minimalSeqIdentity>=10.0) {
            return -2.25;
        }else{
            return -99999.0;
        }
    }else{/*high precision*/
        if (minimalSeqIdentity >=95.0){
            return -0.25;
        } else if (minimalSeqIdentity>=90.0){
            return -0.4;
        } else if (minimalSeqIdentity>=80.0){
            return -0.75;
        }else if (minimalSeqIdentity>=70.0){
            return -1.2;
        }else if (minimalSeqIdentity>=60.0){
            return -1.4;
        }else if (minimalSeqIdentity>=50.0){
            return -1.5;
        }else if (minimalSeqIdentity>=40.0){
            return -1.6;
        }else if (minimalSeqIdentity>=30.0){
            return -1.75;
        } else if (minimalSeqIdentity>=20.0){
            return -1.75;
        }else if (minimalSeqIdentity>=10.0) {
            return -1.8;
        }else{
            return -99999.0;
        }
    }
}/*}}}*/
float GetMinimalSeqIdentity(vector <int> &cntSeqIDTClass, int *maxSeqIDTClass, float* binSeqIDTClass, int numSeqIDTClass){/*{{{*/
    for (int i = 0; i < numSeqIDTClass; i++){
        if (cntSeqIDTClass[i] < maxSeqIDTClass[i]){
            return binSeqIDTClass[2*i];
        }
    }
    return binSeqIDTClass[(numSeqIDTClass-1)*2+1];
}/*}}}*/

struct cmp_str/*{{{*/
{
   bool operator()(char const *a, char const *b)
   {
      return std::strcmp(a, b) < 0;
   }
};/*}}}*/
int CountKMerOfSequence_c(const char* seq, map <char*, int, cmp_str> &wordCountMap, int wordsize){/*{{{*/
/*c style*/
    int len = strlen(seq);
    map <char*,int, cmp_str> :: iterator it;
    int i = 0;
    while (i < len - wordsize + 1 ){
        char *substr=new char[wordsize+1];
        strncpy(substr,seq+i,wordsize);
        substr[wordsize]='\0';
        it = wordCountMap.find(substr);
        if (it == wordCountMap.end()){
            wordCountMap.insert(pair<char*, int > (substr, 1));
        }else {
            it->second += 1;
            delete [] substr;
        }
        i++;
    }
    return wordCountMap.size();
}/*}}}*/
float CalKMerVectorScore_norm_c(map <char*, int, cmp_str> &wordCountMap1, map <char*, int, cmp_str> &wordCountMap2){/*{{{*/
    float bitScore = 0.0;
    float numWord1 = (float)wordCountMap1.size();
    float numWord2 = (float)wordCountMap2.size();
    set <char*, cmp_str> unionKeys;
    map <char*, int, cmp_str> :: iterator it;
    map <char*, int, cmp_str> :: iterator it1;
    map <char*, int, cmp_str> :: iterator it2;
    for (it = wordCountMap1.begin(); it !=  wordCountMap1.end(); it++){
        unionKeys.insert(it->first);
    }
    for (it = wordCountMap2.begin(); it !=  wordCountMap2.end(); it++){
        unionKeys.insert(it->first);
    }
    set <char*, cmp_str> ::iterator iss;
    for (iss = unionKeys.begin(); iss != unionKeys.end(); iss++){
        it1 = wordCountMap1.find(*iss) ;
        it2 = wordCountMap2.find(*iss) ; 
        if (it1 != wordCountMap1.end() && it2 != wordCountMap2.end()){
            bitScore -= fabs(it1->second/numWord1 - it2->second/numWord2);
        }else {
            if (it1 != wordCountMap1.end()){
                bitScore -= it1 -> second/numWord1;
            }else{
                bitScore -= it2 -> second/numWord2;
            }
        }
    }
    return bitScore;
}/*}}}*/
void PrintKMerVectorPair_c(map <char*, int, cmp_str> &wordCountMap1, map <char*, int, cmp_str>& wordCountMap2, const char* id1, const char* id2, FILE *fpout) {/*{{{*/
    set <char*, cmp_str> unionKeys;
    map <char*, int, cmp_str> :: iterator it;
    for (it = wordCountMap1.begin(); it !=  wordCountMap1.end(); it++){
        unionKeys.insert(it->first);
    }
    for (it = wordCountMap2.begin(); it !=  wordCountMap2.end(); it++){
        unionKeys.insert(it->first);
    }
    set <char*, cmp_str> ::iterator iss;
    fprintf(fpout, "\n%13s%10s%13s %d (%d - %d)\n", id1, "", id2, (int)unionKeys.size(), (int)wordCountMap1.size(), (int)wordCountMap2.size());
    for (iss = unionKeys.begin(); iss != unionKeys.end(); iss++){
        it = wordCountMap1.find(*iss) ;
        if (it != wordCountMap1.end()){
            fprintf(fpout,"%6s %6d", (*iss), it->second);
        }else {
            fprintf(fpout,"%6s %6s", (*iss), "NULL");
        }
        fprintf(fpout, "%10s",""); 
        it = wordCountMap2.find(*iss) ;
        if (it != wordCountMap2.end()){
            fprintf(fpout,"%6s %6d", (*iss), it -> second);
        }else {
            fprintf(fpout,"%6s %6s", (*iss), "NULL");
        }
        fprintf(fpout,"\n");
    }
}/*}}}*/
void FreeMem_WrodMap_c(map <char*, int, cmp_str> &wordCountMap){/*{{{*/
    map <char*, int> :: iterator it;
    for (it = wordCountMap.begin(); it !=  wordCountMap.end(); it++){
        delete [] it->first;
    }
}/*}}}*/
float KMerVectorPairwiseComparison_c(string& seq1, string& seq2, string& seqID1, string& seqID2, int wordsize, FILE *fpLog){/*{{{*/
/*c style*/
    map <char*, int, cmp_str> wordCountMap1;
    map <char*, int, cmp_str> wordCountMap2;
    CountKMerOfSequence_c(seq1.c_str(), wordCountMap1, wordsize);
    CountKMerOfSequence_c(seq2.c_str(), wordCountMap2, wordsize);
#ifdef DEBUG_KMER
    PrintKMerVectorPair_c(wordCountMap1, wordCountMap2, seqID1.c_str(), seqID2.c_str(), stdout);
#endif
    float bitScore = CalKMerVectorScore_norm_c(wordCountMap1, wordCountMap2);
    FreeMem_WrodMap_c(wordCountMap1);
    FreeMem_WrodMap_c(wordCountMap2);
    return bitScore;
}/*}}}*/

int CountKMerOfSequence(string &seq, map <string, int> &wordCountMap, int wordsize){/*{{{*/
    int len = seq.size();
    map <string,int> :: iterator it;
    string substr;
    for (int i=0;i < len-wordsize+1; i ++){
        substr = seq.substr(i,wordsize);
        it = wordCountMap.find(substr);
        if (it == wordCountMap.end()){
            wordCountMap.insert(pair<string, int> (substr, 1));
        }else {
            it->second += 1;
        }
    }
    return wordCountMap.size();
}/*}}}*/
void PrintKMerVectorPair(map <string, int> &wordCountMap1, map <string, int>& wordCountMap2, string &id1, string&id2, FILE *fpout) {/*{{{*/
    set <string> unionKeys;
    map <string, int> :: iterator it;
    for (it = wordCountMap1.begin(); it !=  wordCountMap1.end(); it++){
        unionKeys.insert(it->first);
    }
    for (it = wordCountMap2.begin(); it !=  wordCountMap2.end(); it++){
        unionKeys.insert(it->first);
    }
    set <string> ::iterator iss;
    fprintf(fpout, "\n%13s%10s%13s %d (%d - %d)\n", id1.c_str(), "", id2.c_str(), (int)unionKeys.size(), (int)wordCountMap1.size(), (int)wordCountMap2.size());
    for (iss = unionKeys.begin(); iss != unionKeys.end(); iss++){
        it = wordCountMap1.find(*iss) ;
        if (it != wordCountMap1.end()){
            fprintf(fpout,"%6s %6d", (*iss).c_str(), it->second);
        }else {
            fprintf(fpout,"%6s %6s", (*iss).c_str(), "NULL");
        }
        fprintf(fpout, "%10s",""); 
        it = wordCountMap2.find(*iss) ;
        if (it != wordCountMap2.end()){
            fprintf(fpout,"%6s %6d", (*iss).c_str(), it -> second);
        }else {
            fprintf(fpout,"%6s %6s", (*iss).c_str(), "NULL");
        }
        fprintf(fpout,"\n");
    }
}/*}}}*/
float CalKMerVectorScore_old1(map <string, int> &wordCountMap1, map <string, int> &wordCountMap2){/*{{{*/
    float bitScore = 0.0;
    set <string> unionKeys;
    map <string, int> :: iterator it;
    map <string, int> :: iterator it1;
    map <string, int> :: iterator it2;
    for (it = wordCountMap1.begin(); it !=  wordCountMap1.end(); it++){
        unionKeys.insert(it->first);
    }
    for (it = wordCountMap2.begin(); it !=  wordCountMap2.end(); it++){
        unionKeys.insert(it->first);
    }
    set <string> ::iterator iss;
    for (iss = unionKeys.begin(); iss != unionKeys.end(); iss++){
        it1 = wordCountMap1.find(*iss) ;
        it2 = wordCountMap2.find(*iss) ; 
        if (it1 != wordCountMap1.end() && it2 != wordCountMap2.end()){
            bitScore += min(it1->second, it2->second);
        }else {
            if (it1 != wordCountMap1.end()){
                bitScore -= it1 -> second;
            }else{
                bitScore -= it2 -> second;
            }
        }
    }
    return bitScore;
}/*}}}*/
float CalKMerVectorScore(map <string, int> &wordCountMap1, map <string, int> &wordCountMap2){/*{{{*/
    float bitScore = 0.0;
    set <string> unionKeys;
    map <string, int> :: iterator it;
    map <string, int> :: iterator it1;
    map <string, int> :: iterator it2;
    for (it = wordCountMap1.begin(); it !=  wordCountMap1.end(); it++){
        unionKeys.insert(it->first);
    }
    for (it = wordCountMap2.begin(); it !=  wordCountMap2.end(); it++){
        unionKeys.insert(it->first);
    }
    set <string> ::iterator iss;
    for (iss = unionKeys.begin(); iss != unionKeys.end(); iss++){
        it1 = wordCountMap1.find(*iss) ;
        it2 = wordCountMap2.find(*iss) ; 
        if (it1 != wordCountMap1.end() && it2 != wordCountMap2.end()){
            bitScore -= abs(it1->second - it2->second);
        }else {
            if (it1 != wordCountMap1.end()){
                bitScore -= it1 -> second;
            }else{
                bitScore -= it2 -> second;
            }
        }
    }
    return bitScore;
}/*}}}*/
float CalKMerVectorScore_norm(map <string, int> &wordCountMap1, map <string, int> &wordCountMap2){/*{{{*/
    float bitScore = 0.0;
    float numWord1 = (float)wordCountMap1.size();
    float numWord2 = (float)wordCountMap2.size();
    set <string> unionKeys;
    map <string, int> :: iterator it;
    map <string, int> :: iterator it1;
    map <string, int> :: iterator it2;
    for (it = wordCountMap1.begin(); it !=  wordCountMap1.end(); it++){
        unionKeys.insert(it->first);
    }
    for (it = wordCountMap2.begin(); it !=  wordCountMap2.end(); it++){
        unionKeys.insert(it->first);
    }
    set <string> ::iterator iss;
    for (iss = unionKeys.begin(); iss != unionKeys.end(); iss++){
        it1 = wordCountMap1.find(*iss) ;
        it2 = wordCountMap2.find(*iss) ; 
        if (it1 != wordCountMap1.end() && it2 != wordCountMap2.end()){
            bitScore -= fabs(it1->second/numWord1 - it2->second/numWord2);
        }else {
            if (it1 != wordCountMap1.end()){
                bitScore -= it1 -> second/numWord1;
            }else{
                bitScore -= it2 -> second/numWord2;
            }
        }
    }
    return bitScore;
}/*}}}*/
float KMerVectorPairwiseComparison(string& seq1, string& seq2, string& seqID1, string& seqID2, int wordsize, FILE *fpLog){/*{{{*/
    map <string, int> wordCountMap1;
    map <string, int> wordCountMap2;
    CountKMerOfSequence(seq1, wordCountMap1, wordsize);
    CountKMerOfSequence(seq2, wordCountMap2, wordsize);
#ifdef DEBUG_KMER
    PrintKMerVectorPair(wordCountMap1, wordCountMap2, seqID1, seqID2, stdout);
#endif
    //return 0.0;
    //float bitScore = CalKMerVectorScore(wordCountMap1, wordCountMap2);
    float bitScore = CalKMerVectorScore_norm(wordCountMap1, wordCountMap2);
    //bitScore = bitScore*(max(len1,len2)/(min(len1,len2) + 1e-9));
    return bitScore;
}/*}}}*/


void PrintKMerVectorPair_table(int* freqTable1, int* freqTable2,  vector<int>& indexWordHit1, vector <int> & indexWordHit2, string& id1, string& id2, FILE *fpout) {/*{{{*/
    set <int> unionKeys;
    unionKeys.insert(indexWordHit1.begin(),indexWordHit1.end());
    unionKeys.insert(indexWordHit2.begin(),indexWordHit2.end());
    set <int> ::iterator it;
    fprintf(fpout, "\n%13s%10s%13s %d (%d - %d)\n", id1.c_str(), "", id2.c_str(), (int)unionKeys.size(), (int)indexWordHit1.size(), (int)indexWordHit2.size());
}/*}}}*/
int CountKMerOfSequence_table(const char*seq, int len, int* freqTable, set <int>&indexWordHit, int wordsize){/*{{{*/
    set <int> ::iterator it;
    for (int i= 0; i < len - wordsize + 1 ; i++){
        int idx = 0;
        for (int j = 0; j < wordsize; j++){
            idx +=  *(aacode+seq[i+j]) * para20[wordsize-j-1];
        }
        it = indexWordHit.find(idx);
        if (it == indexWordHit.end()){
            indexWordHit.insert(idx);
            freqTable[idx] = 1;
        }else{
            freqTable[idx] ++;
        }
    }
    return indexWordHit.size();
}/*}}}*/
int CountKMerOfSequence_table_2(const char*seq, int len, int *freqTable, vector <int>&indexWordHit, int wordsize){/*{{{*/
    int numTotalWord = len-wordsize +1; /*number of words with duplicates*/
    int *idxList = new int [numTotalWord];
    for (int i= 0; i < numTotalWord ; i++){
        int idx = 0;
        for (int j = 0; j < wordsize; j++){
            idx +=  *(aacode+seq[i+j]) * para20[wordsize-j-1];
        }
        idxList[i] = idx;
        freqTable[idx] = 0;
    }
    for (int i= 0; i < numTotalWord ; i++){
        freqTable[idxList[i]] ++;
    }
    vector <int> newlist(idxList, idxList+numTotalWord);
    vector <int> :: iterator it;
    sort (newlist.begin(),newlist.end());
    it = unique(newlist.begin(),newlist.end()); 
    newlist.resize(it - newlist.begin());/*this must be added to actually get a uniq list*/
    indexWordHit = newlist;
    //indexWordHit.insert(idxList, idxList+numTotalWord); [>this takes 90% of time for this function<]
    delete [] idxList;
    return indexWordHit.size();
}/*}}}*/
float CalKMerVectorScore_norm_table(int * freqTable1, int* freqTable2, set <int> indexWordHit1, set <int> indexWordHit2){/*{{{*/
    float bitScore = 0.0;
    float numWord1 = (float)indexWordHit1.size();
    float numWord2 = (float)indexWordHit2.size();
    set <int> index_intersection;
    set <int> index_1_minus_2;
    set <int> index_2_minus_1;
    set_intersection(indexWordHit1.begin(), indexWordHit1.end(),
            indexWordHit2.begin(),indexWordHit2.end(),
            inserter(index_intersection,index_intersection.begin()));
    set_difference(indexWordHit1.begin(), indexWordHit1.end(),
            indexWordHit2.begin(),indexWordHit2.end(), 
            inserter(index_1_minus_2, index_1_minus_2.begin()));
    set_difference(indexWordHit2.begin(), indexWordHit2.end(),
            indexWordHit1.begin(),indexWordHit1.end(), 
            inserter(index_2_minus_1, index_2_minus_1.begin()));

    set <int> ::iterator it;
    for (it = index_intersection.begin(); it != index_intersection.end(); it++){
        bitScore -= fabs(freqTable1[*it]/numWord1 - freqTable2[*it]/numWord2);
    }
    for (it = index_1_minus_2.begin(); it != index_1_minus_2.end(); it++){
        bitScore -= fabs(freqTable1[*it]/numWord1);
    }
    for (it = index_2_minus_1.begin(); it != index_2_minus_1.end(); it++){
        bitScore -= fabs(freqTable2[*it]/numWord2);
    }
    return bitScore;
}/*}}}*/
float CalKMerVectorScore_norm_table_2(int * freqTable1, int* freqTable2, vector <int> indexWordHit1, vector <int> indexWordHit2){/*{{{*/
    float bitScore = 0.0;
    float numWord1 = (float)indexWordHit1.size();
    float numWord2 = (float)indexWordHit2.size();
    vector <int> index_intersection;
    vector <int> index_1_minus_2;
    vector <int> index_2_minus_1;
    set_intersection(indexWordHit1.begin(), indexWordHit1.end(),
            indexWordHit2.begin(),indexWordHit2.end(),
            inserter(index_intersection,index_intersection.begin()));
    set_difference(indexWordHit1.begin(), indexWordHit1.end(),
            indexWordHit2.begin(),indexWordHit2.end(), 
            inserter(index_1_minus_2, index_1_minus_2.begin()));
    set_difference(indexWordHit2.begin(), indexWordHit2.end(),
            indexWordHit1.begin(),indexWordHit1.end(), 
            inserter(index_2_minus_1, index_2_minus_1.begin()));

#ifdef DEBUG_SET_OPT/*{{{*/
    fprintf(stdout,"vector1(%d):\n", int(numWord1));
    for (int i = 0; i < int(numWord1); i++){
        fprintf(stdout,"%3d ", indexWordHit1[i]);
    }
    fprintf(stdout,"\n");

    fprintf(stdout,"vector2(%d):\n", int(numWord2));
    for (int i = 0; i < int(numWord2); i++){
        fprintf(stdout,"%3d ", indexWordHit2[i]);
    }
    fprintf(stdout,"\n");

    fprintf(stdout,"vector intersection(%d):\n", int(index_intersection.size()));
    for (int i = 0; i < index_intersection.size(); i++){
        fprintf(stdout,"%3d ", (index_intersection.begin())[i]);
    }
    fprintf(stdout,"\n");

    fprintf(stdout,"vector diff 1(%d):\n", int(index_1_minus_2.size()));
    for (int i = 0; i < index_1_minus_2.size(); i++){
        fprintf(stdout,"%3d ", (index_1_minus_2.begin())[i]);
    }
    fprintf(stdout,"\n");
#endif /*}}}*/

    vector <int> ::iterator it;
    for (it = index_intersection.begin(); it != index_intersection.end(); it++){
        bitScore -= fabs(freqTable1[*it]/numWord1 - freqTable2[*it]/numWord2);
    }
    for (it = index_1_minus_2.begin(); it != index_1_minus_2.end(); it++){
        bitScore -= (freqTable1[*it]/numWord1);
    }
    for (it = index_2_minus_1.begin(); it != index_2_minus_1.end(); it++){
        bitScore -= (freqTable2[*it]/numWord2);
    }
    return bitScore;
}/*}}}*/
float KMerVectorPairwiseComparison_table(string& seq1, string& seq2, string& seqID1, string& seqID2, int* freqTable1, int *freqTable2, int sizeTable, int wordsize, FILE *fpLog){/*{{{*/
    /*keep a table*/
    int len1 = seq1.size();
    int len2 = seq2.size();
    vector <int> indexWordHit1;
    vector <int> indexWordHit2;
    CountKMerOfSequence_table_2(seq1.c_str(), len1, freqTable1, indexWordHit1, wordsize);
    CountKMerOfSequence_table_2(seq2.c_str(), len2, freqTable2, indexWordHit2, wordsize);
#ifdef DEBUG_KMER
    PrintKMerVectorPair_table(freqTable1,freqTable2, indexWordHit1,indexWordHit2, seqID1, seqID2, stdout);
#endif
    //return 0.0;
    float bitScore ; 
    //float bitScore = CalKMerVectorScore_norm_table(freqTable1, freqTable2, indexWordHit1, indexWordHit2);
    bitScore = CalKMerVectorScore_norm_table_2(freqTable1, freqTable2, indexWordHit1, indexWordHit2);
    //bitScore = bitScore*(max(len1,len2))/1000;
    return bitScore;
}/*}}}*/

int ReadInPairList(const char* file, vector <Pair>&pairList)/*{{{*/
{
    ifstream ifp (file, ios::binary);
    if (ifp.is_open()){
        ifp.seekg (0, ios::end);
        int length = ifp.tellg();
        ifp.seekg (0, ios::beg);
        char *buffer = new char [ length+1];
        ifp.read(buffer,length);
        ifp.close();
        buffer[length] = '\0';
        char *pch;
        pch = strtok(buffer, "\n");
        Pair tmppair;
        int status_sscanf;
        while (pch != NULL){
            int linesize=strlen(pch);
            if ( linesize > 0 ){
                char *mem1=new char[linesize];
                char *mem2=new char[linesize];
                if ((status_sscanf = sscanf(pch, "%s %s", mem1, mem2)) == 2){
                    tmppair.mem1 = mem1;
                    tmppair.mem2 = mem2;
                    pairList.push_back(tmppair);
                }
                delete [] mem1;
                delete [] mem2;
            }
            pch = strtok(NULL, "\n");
        }
        delete [] buffer;
        return pairList.size();
    }else{
        cerr << "Failed to open pairfile "
            << file 
            << endl;
        return -1;
    }
    return 0;

}/*}}}*/
char *SpanExcluding(const char* strForSpan,char* strAfterSpan, const char charSet[]/*=WHITE_SPACE*/)/*{{{*/
//**********************************************************************
// SpanExcluding()
// compress the string by removing the characters in the charSet
// and store the compressed string into strAfterSpan
// the search is case-sensitive
// default charSet is WHITE_SPACE
//**********************************************************************
{
	int l = strlen(strForSpan);
	int n = strlen(charSet);
	int i,j ; 
	j = 0 ;
	for(i = 0 ; i < l ;i ++)
	{
		if(!IsInCharSet(strForSpan[i],charSet, n))
		{
			strAfterSpan[j] = strForSpan[i] ;
			j ++ ;
		}
	}
	strAfterSpan[j] = '\0' ;
	return strAfterSpan;
}/*}}}*/
int my_strcpy(char* to, const char* from, int max, int sizefrom = 0)/*{{{*/
/******************************************************************************
 * my modification of strcpy
 * copy max characters from "from" to "to", add NULL terminator automatically
 * updated 2008-04-23, memcpy win in speed when the copying string is long
 * updated 2011-10-27:
 *   since calling strlen will be time consuming for very large from string,
 *   e.g. when "from" is the buffer of a whole trunk of file. Therefore, strlen
 *   part is removed.
 *****************************************************************************/
{
    if(sizefrom < 200) {
        strncpy(to,from,max);
    } else {
        memcpy(to,from,max);
    }
    to[max] = '\0';
    return max;
}/*}}}*/
int checkfilestream(FILE *fp, const char* filename, const char *mode, bool isAssert =false)/*{{{*/
{
    if( fp == NULL) {
        fprintf(stderr,"Can not open file '%s' with mode '%s'\n", filename,mode);
        if(isAssert) {
            assert(fp != NULL);
        }
        return -1;
    }else{
        return 0;
    }
}
/*}}}*/
char *rootname(const char* filename, char* rtname, int max_rtname )/*{{{*/
/*****************************************************************************
 * rootname
 * given the file name, 
 * return the rootname of the filename
 ****************************************************************************/
{
    const char *pch;
    char *pstr;
    if((pch = strrchr(filename,'/')) != NULL)
        pstr = (char*) pch+1;
    else
        pstr = (char*) filename;

    if((pch = strrchr(pstr,'.')) != NULL)
        my_strcpy(rtname,pstr, min((int)(pch - pstr), max_rtname));
    else
        rtname = pstr;
    return rtname;
}
/*}}}*/
int fgetline(FILE* fp, char* line, int max/* = 0x7FFFFFFF*/)/*{{{*/
/*****************************************************************************
 * Read one line from fp, copying it to line array (but no more than max
 * chars). Does not place terminating \n in line array.  
 * it can be called without "max" flag, but should make sure the allocated
 * memory for "line" is larger than the longest line.
 *
 * Returns: line length, or 0 for empty line, or "EOF" for end-of-file.
 * 
 * since fgetc(), getc(), getchar(), returns int value,always use an int
 * variable to store the result of the fgetc(), getc() and getchar().  
 * getc() is faster than fgetc()
 *   LOG: 2006-04-26 16:31:29 Wednesday  Week 17 <nanjiang@shu>
 *   bug fixed for '\n' return keys, 
 *   now it is valid for both dos and unix
 ****************************************************************************/
{
    int nch = 0; /* record number of characters actually read */
    int c;
    max = max - 1;			/* leave room for '\0' */
    while((c = getc(fp)) != EOF)
    {
        if(c == 0x0d ) continue; /* in unix, '\n'= 0x0a, thus 0x0d will be cheated as another character*/ 
        if(c == '\n') break; /* in dos, '\n' is also equal to 0x0a, but the preceding 0x0d will not be read*/ 

        if(nch < max)
        {
            line[nch] = c;
            nch = nch + 1;
        }
    }
    line[nch] = '\0';

    if(c == EOF && nch == 0) return EOF;
    else return nch;
}/*}}}*/
static void printPathMatrix(const float* path, const ajint* compass,/*{{{*/
	const char *a, const char *b, ajuint lena, ajuint lenb, FILE *fpLog)
{
    char compasschar;
    ajuint i;
    ajuint j;

    fprintf(fpLog,"path matrix:\n");

    i = lena;
    while( i--!= 0)
    {
        fprintf(fpLog, "%4d(%c)", i, a[i]);

        for(j = 0; j < lenb; j++)
        {
            if(compass[i * lenb + j] == LEFT)
                compasschar = '<';
            else if(compass[i * lenb + j] == DOWN)
                compasschar = '^';
            else
                compasschar = ' ';

            fprintf(fpLog, "%6.2f%c ", path[i * lenb + j], compasschar);
        }
        fprintf(fpLog,"\n");
    }

    fprintf(fpLog,"       ");

    for (j = 0; j < lenb; ++j)
        fprintf(fpLog, "%4d(%c) ", j, b[j]);

    fprintf(fpLog, "\n");

    return;
}/*}}}*/
int Char2Digit(char aa, const char* alphabet, int n = 0)/*{{{*/
/****************************************************************************
 * Char2Digit()
 *
 * return the digital indexing of the character in alphabet, 
 * starting from 0
 * if the character is not in the alphabet, return -1
 ***************************************************************************/
{
    if( n == 0) 
        n = strlen(alphabet);
    const char *pch;
    if((pch = strchr(alphabet,aa))!= NULL)
    {
        int index = pch - alphabet;
        return index;
    } else {
        return (n-1); // if not in alphabet, consider as * or unknown, which is the last char in the alphabet
    }
}/*}}}*/
int Charcase2Digit(char aa, const char* alphabet, int n /*= 0*/)/*{{{*/
/****************************************************************************
 * Charcase2Digit()
 *
 * return the digital indexing of the character in alphabet, case insensitive, 
 * starting from 0
 * if the character is not in the alphabet, return -1
 ***************************************************************************/
{
    if( n == 0) {
        n = strlen(alphabet);
    }
    const char *pch;
    char aa_lower = tolower(aa);
    char aa_upper = toupper(aa);
    if((pch = strchr(alphabet,aa_lower))!= NULL 
            ||(pch = strchr(alphabet, aa_upper))!= NULL) {
        int index = pch - alphabet;
        return index;
    } else {
        return (n-1); // if not in alphabet
    }
}/*}}}*/
char * StrReverse(char *str)/*{{{*/
{
    if (!str || strcmp(str,"")== 0)
    {
        return str;
    }
    else 
    {
        int n = strlen(str);
        Array1D <char> tmpstr(n+10);
        strcpy(tmpstr.array1D,str);
        int i;
        for (i = 0; i < n ; i++)
        {
            str[i] = tmpstr.array1D[n-i-1];
        }
        return str;
    }
}/*}}}*/
template <class T> int ReadSMatrix(const char *filename, T **S, char *alphabet)/*{{{*/
/*****************************************************************************
 * read in substitution matrix for protein sequence alignment, 
 * e.g., BLOSUM62, which is also the default substitution matrix
 ****************************************************************************/
{
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 

    char *pch;
    int i,j;
    FILE* fpin = NULL;
    fpin = fopen(filename,"r");
    if (fpin == NULL){
        fprintf(stderr,"Failed to read file %s\n", filename);
        return -1;
    }else{
        while(fgetline(fpin, line,maxline) != EOF) {
            if(line[0] != '#') break ;
        }
        SpanExcluding(line,alphabet," ");
        int n = strlen(alphabet);
        char delim[] = WHITE_SPACE;
        i = 0;
        while((linesize = fgetline(fpin, line,maxline))!= EOF) {
            pch = strtok(line, delim);
            j = 0;
            while(pch != NULL) {
                if(j >= 1 && j < n+1) S[i][j-1] = T (atof(pch));
                j ++;
                pch = strtok(NULL,delim);
            }
            i ++;
        }
        fclose(fpin);
        if (i != n) {
            fprintf(stderr,"Error! The number of rows and columns in the matrix %s are not equal.\n",filename);
            return -1;
        }
        return n;
    }
}
template int ReadSMatrix<int>   (const char* filename, int **S   , char *alphabet);
template int ReadSMatrix<float> (const char* filename, float **S , char *alphabet);
template int ReadSMatrix<double>(const char* filename, double **S, char *alphabet);
/*}}}*/
template <class T> void Array1Dto2D(T *array1d, T ** array2d, int m, int n)/*{{{*/
{
    int i,j,idx;
    for (i = 0; i <m; i++ )
    {
        for (j = 0; j <n; j++ )
        {
            idx = i * n + j;
            array2d[i][j] = array1d[idx];
        }
    }
}/*}}}*/
int ReadNextSeq_FASTA(FILE *fp, char* seq, float *gapopenArray, float *gpeArray, float *tgpeArray, bool &isHasGapOpenArray, int *pSeq_type /*= NULL*/, int maxlength /*= LONGEST_SEQ*/, char* annotationLine/*=NULL*/, int maxSizeAnnotationLine /*=50*/)/*{{{*/
    /****************************************************************************
     * ReadNextSeq_FASTA()
     * read in the fasta format sequence from the file stream 
     * return the length of the sequence if successful
     * return 0 or minus value if no more sequence can be read from the file stream
     * The leading white spaces are ignored
     * check the type (DNA or AA) of fasta sequence file based the annotation line 
     * Last modified, 2007-02-12, Nanjiang Shu
     * Updated 2010-09-09
     * The gapopen values are read in if supplied, the gapopen values are
     * enclosed in {  }
     *
     * Updated 2010-10-07 
     * #1. tags gpo gapopen
     *       gpe gapextension
     *       tgpe terminal gapextension 
     *  can be supplied, if tag is not supplied, assumed to be gpo
     *  2. { #xxxxx  } the leading '#' means this gap penalty array is commented out
     ***************************************************************************/
{
    int   c;
    int   i;
    int cntVar = 0;
    isHasGapOpenArray = false;
    do{  /* ignore the leading white spaces */ 
        c = getc(fp);
    }while(isspace(c));

    if(pSeq_type != NULL) 
        *pSeq_type = UNKNOWN_SEQ_TYPE; /* initializing sequence type*/

    if(c  == '>') 
    { 
        Array1D <char> line_1darray(maxSizeAnnotationLine +1);
        char *line = line_1darray.array1D;
        fgetline(fp,line,maxSizeAnnotationLine);
        if( pSeq_type != NULL)
        {
            if(strstr(line, "protein") != NULL)
                *pSeq_type = AA_SEQ;
            else if( strstr(line,"nucleic") != NULL)
                *pSeq_type = DNA_SEQ;
            else if( strstr(line,"shape") != NULL)
                *pSeq_type = SHAPE_SEQ;
            else
                *pSeq_type = UNKNOWN_SEQ_TYPE;
        }
        if (annotationLine != NULL)
        {
            my_strcpy(annotationLine,line, maxSizeAnnotationLine);  /*read in the annotation line*/ 
        }
    }
    else  
    {
        fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream if there is no annotation line*/ 
    }

    i = 0 ;
    while((c = getc(fp)) != EOF) 
    {
        if(c == '>') 
        {
            fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream*/
            break;  /* read in the first sequence if there are multiple sequences in the file*/ 
        }else if (c == '{'){
            while(isspace(c = getc(fp))); /*neglect leading spaces*/  
            if (c == '#'){ /*if the leadning nonblank char is #, neglect*/  
                while((c=getc(fp)) != '}');
            }else if (c != '}' ){
                char tmpstr[100] = "";
                int cntdigit = 0;
                float *pg = NULL;
                cntVar = 0;
                isHasGapOpenArray = true;
                cntdigit = 0;
                fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream*/
                while(!isspace(c=getc(fp))) { /*get the first word*/
                    tmpstr[cntdigit] = c;
                    cntdigit++;
					if (c == ':'){
						break;
					}
                    if (cntdigit > 100) {
                        fprintf(stderr,"word too long in gapopen arry in fasta file. WORD=%100s\ncode line %d file %s\n", tmpstr, __LINE__, __FILE__);
                        exit(1);
                    }
                }
                tmpstr[cntdigit] = '\0';
#ifdef DEBUG
                printf("function %s: The first word is:%s\n", __FUNCTION__, tmpstr);
#endif 
                if (strcmp(tmpstr,"gpo:") == 0 || strcmp(tmpstr,"gapopen:")==0){
                    pg = gapopenArray;
                }else if (strcmp(tmpstr,"gpe:") == 0){
                    pg = gpeArray;
                }else if (strcmp(tmpstr,"tgpe:") == 0){
                    pg = tgpeArray;
                }else if (isdigit((int)tmpstr[0])){ /*if no tag is supplied, by default it is gpo*/
                    pg = gapopenArray;
                    if (sscanf(tmpstr,"%f", &pg[cntVar]) == 1){
                        cntVar++;
                    }
                }
                cntdigit=0;
                while((c = getc(fp)) != '}'){
                    if (!isspace(c)){
                        tmpstr[cntdigit] = c;
                        cntdigit ++;
                    } else if (cntdigit>0){
                        tmpstr[cntdigit] = '\0';
                        if (sscanf(tmpstr,"%f", &pg[cntVar]) == 1){
                            cntVar++;
                        }
                        cntdigit = 0;
                    }
                }
				tmpstr[cntdigit] = '\0';
                if (strcmp(tmpstr, "") != 0){
                    if (sscanf(tmpstr,"%f", &pg[cntVar]) == 1){
                        cntVar++;
                    }
                }
            }
        }else if (isalpha(c)){ /* neglect white spaces and return characters in sequence region*/
            seq[i] = c ; i ++ ;
            if(i >= maxlength){
                fprintf(stderr,"Error, sequence longer then maxlength = %d\n", maxlength);
                exit(1);
            }
        }
    }
    seq[i] = '\0' ;
    if (isHasGapOpenArray && i != cntVar)
    {
        fprintf(stderr,"Error, the seqLength (%d) and the number of gapopen values (%d) do not match\n", i, cntVar);
        return -1;
    }
    if(c == EOF && i == 0)
        return EOF;
    else
        return i; /* return the length of sequence*/ 
}/*}}}*/
int ReadSeq_FASTA(const char *fileName, char* seq, float *gapopenArray, float * gpeArray, float*tgpeArray, bool &isHasGapOpenArray, int *pSeq_type /*= NULL*/, int maxlength /*= LONGEST_SEQ*/, char* annotationLine/*=NULL*/, int maxSizeAnnotationLine /*=50*/)/*{{{*/
    /****************************************************************************
     * ReadSeq_FASTA()
     * read in fasta format sequence file which contains single sequence
     * return the length of the sequence if successful
     * return -1 if can not open file

     * LOG: 2006-04-26 16:27:02 Wednesday  Week 17 <nanjiang@shu>
     * bug fixed for 0x0d character in unix system, in windows, '\n' is 0x0d0x0a,
     * while under unix, '\n' = 0x0a;
     * also the leading white spaces are neglected
     * LOG: 2006-06-14 16:36:56 Wednesday  Week 24 <nanjiang@casio>
     *   check the type (DNA or AA) of fasta sequence file, using the annotation line 
     * Updated 2010-04-22: the annotation line (without the leading ">") can be
     *   read in, note that the maxSizeAnnotationLine must be supplied
     * Updated 2010-09-09
     * The gapopen values are read in if supplied, the gapopen values are
     * enclosed in {  }
     * Updated 2010-10-07 
     * #1. tags gpo gapopen
     *       gpe gapextension
     *       tgpe terminal gapextension 
     *  can be supplied, if tag is not supplied, assumed to be gpo
     *  2. { #xxxxx  } the leading '#' means this gap penalty array is commented out
     ***************************************************************************/
{
    int   c;
    int   i;
    FILE *fp = NULL;
    int cntVar = 0;
    isHasGapOpenArray = false;
    if((fp = fopen(fileName,"r")) == NULL){
        fprintf(stderr, "can not open file %s \n",fileName) ;
        exit(1);
    }

    do{  /* neglect the leading white spaces */ 
        c = getc(fp);
    }while(isspace(c));

    if(pSeq_type != NULL) 
        *pSeq_type = UNKNOWN_SEQ_TYPE; /* initializing sequence type*/

    if(c  == '>'){ 
        Array1D <char> line_1darray(maxSizeAnnotationLine +1);
        char *line = line_1darray.array1D;
        fgetline(fp,line,maxSizeAnnotationLine);
        if( pSeq_type != NULL)
        {
            if(strstr(line, "protein") != NULL)
                *pSeq_type = AA_SEQ;
            else if( strstr(line,"nucleic") != NULL)
                *pSeq_type = DNA_SEQ;
            else if( strstr(line,"shape") != NULL)
                *pSeq_type = SHAPE_SEQ;
            else
                *pSeq_type = UNKNOWN_SEQ_TYPE;
        }
        if (annotationLine != NULL) {
            my_strcpy(annotationLine,line, maxSizeAnnotationLine);  /*read in the annotation line*/ 
        }
    }else{
        fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream if there is no annotation line*/ 
    }

    i = 0 ;
    while((c = getc(fp)) != EOF){
        if(c == '>'){
            break;  /* read in the first sequence if there are multiple sequences in the file*/ 
        }else if (c == '{'){
            while(isspace(c = getc(fp))); /*neglect leading spaces*/  
            if (c == '#'){ /*if the leadning nonblank char is #, neglect*/  
                while((c=getc(fp)) != '}');
            } else if (c != '}'){
                int cntdigit = 0;
                char tmpstr[100] = "";
                float *pg = NULL;
                cntVar = 0;
                isHasGapOpenArray = true;
                fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream*/
                while(!isspace(c=getc(fp))) { /*get the first word*/
                    tmpstr[cntdigit] = c;
                    cntdigit++;
					if (c == ':'){
						break;
					}
                    if (cntdigit > 100) {
                        fprintf(stderr,"word too long in gapopen arry in fasta file. WORD=%100s\ncode line %d file %s\n", tmpstr, __LINE__, __FILE__);
                        exit(1);
                    }
                }
                tmpstr[cntdigit] = '\0';
#ifdef DEBUG
                printf("function:%s, first word=%s\n", __FUNCTION__, tmpstr);
#endif
                if (strcmp(tmpstr,"gpo:") == 0 || strcmp(tmpstr,"gapopen:")==0){
                    pg = gapopenArray;
                }else if (strcmp(tmpstr,"gpe:") == 0){
                    pg = gpeArray;
                }else if (strcmp(tmpstr,"tgpe:") == 0){
                    pg = tgpeArray;
                }else if (isdigit((int)tmpstr[0])){ /*if no tag is supplied, by default it is gpo*/
                    pg = gapopenArray;
                    if (sscanf(tmpstr,"%f", &pg[cntVar]) == 1){
                        cntVar++;
                    }
                }
                else{
                    fprintf(stderr,"wrong gap penalty tag in file %s:\"%s\"\n", fileName, tmpstr);
                    exit(1);
                }
                cntdigit=0;
                while((c = getc(fp)) != '}'){
                    if (!isspace(c)){
                        tmpstr[cntdigit] = c;
                        cntdigit ++;
                    } else if (cntdigit>0){
                        tmpstr[cntdigit] = '\0';
                        if (sscanf(tmpstr,"%f", &pg[cntVar]) == 1){
                            cntVar++;
                        }
                        cntdigit = 0;
                    }
                }
				tmpstr[cntdigit] = '\0';
                if (strcmp(tmpstr, "") != 0){
                    if (sscanf(tmpstr,"%f", &pg[cntVar]) == 1){
                        cntVar++;
                    }
                }
            }
        }else if (isalpha(c)){ /* neglect white spaces and return characters in sequence region*/
            seq[i] = c ; i ++ ;
            if(i >= maxlength){
                fprintf(stderr,"Error, sequence longer then maxlength = %d\n", maxlength);
                exit(1);
            }
        }
    }
    seq[i] = '\0' ;
    if (isHasGapOpenArray && i != cntVar)
    {
        fprintf(stderr,"Error, the seqLength (%d) and the number of gapopen values (%d) do not match\n", i, cntVar);
        return -1;
    }
    fclose(fp);
    return i; /* return the length of sequence*/ 
}/*}}}*/
string GetIDFromAnnotationLine(char *line, string &id)/*{{{*/
{
    char tmpstr[100] ="";
    if (line[0] == '>') {
        line[0] = ' ';
    }
    sscanf(line, " %100s ", tmpstr );
    char * pchp;
    char * pch;
    pchp = &tmpstr[0]-1;
    pch = strchr(tmpstr,'|');
    while (pch != NULL) {
        pchp = pch;
        pch=strchr(pch+1,'|');
    }
    id = pchp+1;
    return id;
}/*}}}*/
int ReadNextGapOpenArray(FILE *fp, float *gapopenArray)/*{{{*/
/*read in file in the format
 * { 10 10 10 }
 * { 20 20 20 }
 * */
{
    int cntVar = 0;
    int c = 0;
    while((c = getc(fp)) != EOF){
        if (c == '{') { 
            char tmpstr[100] = "";
            int cntdigit = 0;
            while((c = getc(fp)) != '}') {
                if (!isspace(c)) {
                    tmpstr[cntdigit] = c;
                    cntdigit ++;
                } else {
                    tmpstr[cntdigit] = '\0';
                    if (sscanf(tmpstr,"%f", &gapopenArray[cntVar]) == 1) {
                        cntVar++;
                    }
                    strcpy(tmpstr,"");
                    cntdigit = 0;
                }
            }
            if (strcmp(tmpstr, "") != 0) {
                sscanf(tmpstr,"%f", &gapopenArray[cntVar]);
                cntVar++;
            }
            break;
        }       
    }
    return cntVar;
}/*}}}*/
int GetSubstitutionMatrix(string &file, float **SM, char * alphabet, int alnType)/*{{{*/
{
    if (file != "" ) {
        if (ReadSMatrix(file.c_str(), SM, alphabet) == -1){
            return -1;
        }
    } else {
        int i,j;
        if (alnType == 0) {
            for (i=0;i<NUM_NUC;i++){
                for(j=0;j<NUM_NUC;j++){
                    SM[i][j] = float(nuc44[i][j]);
                }
            }
            strcpy(alphabet, ALPHABET_NUC);
        } else if (alnType == 1) {
            for (i=0;i<NUM_BLOSUM;i++){
                for(j=0;j<NUM_BLOSUM;j++){
                    SM[i][j] = float(blosum62[i][j]);
                }
            }
            strcpy(alphabet, ALPHABET_BLOSUM);
        } else if (alnType == 2) {
            for (i=0;i<NUM_SHAG;i++){
                for(j=0;j<NUM_SHAG;j++){
                    SM[i][j] = float(shag6[i][j]);
                }
            }
            strcpy(alphabet, ALPHABET_SHAG);
        } else {
            fprintf(stderr,"Wrong alignment type:'%d'\n", alnType);
            return -1;
        }
    }
    return 0;
}/*}}}*/
void AlignAna(const char *alignX,const char* alignY, const char *alphabet, int length, float **subMatr, AlignFactor* pAlignFactor, int* alignRel)/*{{{*/
{
    int i ;
    int matchCnt = 0;
    int simiCnt  = 0;
    int gapCnt   = 0;
    for(i = length-1 ; i >= 0 ; i -- )
    {
        if(alignX[i] == CHAR_GAP || alignY[i] == CHAR_GAP)
        {
            gapCnt ++ ;
            alignRel[i] = GAP ;
        }
        else
        {
            if(alignX[i] == alignY[i])
            {
                matchCnt ++ ;
                alignRel[i] = IDT ;
            }
            else if(subMatr[Char2Digit(alignX[i],alphabet)][Char2Digit(alignY[i], alphabet)] > 0)
            {
                simiCnt ++ ;
                alignRel[i] = SIM ;
            }
            else
            {
                alignRel[i] = MIS ;
            }
        }
    }
    simiCnt += matchCnt;
    pAlignFactor->idt_cnt    = matchCnt;
    pAlignFactor->sim_cnt    = simiCnt;
    pAlignFactor->gap_cnt    = gapCnt;
    pAlignFactor->identity   = float(matchCnt) / length;
    pAlignFactor->similarity = float((simiCnt - matchCnt)*1.0 + matchCnt) / length;
    pAlignFactor->gapPercent = float(gapCnt)   / length;
}
/*}}}*/
void DebugPrintPathMatrix2(const char *header, float ** V, char *Xstr,char *Ystr,int m,int n,int **compass, FILE *fplog)/*{{{*/
{
    fprintf(fplog, "%s\n", header);
    int i,j;
    for( i = m-1 ; i >= 0 ; i --)
    {   
        fprintf(fplog,"%4d(%c)",i, Xstr[i]);
        for(j= 0 ; j < n ; j ++)
        {
            char tr = ' ';
            if (compass[i][j] == DOWN)
            {
                tr = '^';
            }
            else if (compass[i][j] == LEFT)
            {
                tr = '<';
            }
            fprintf(fplog,"%6.2f%c ",V[i][j], tr);
        }
        fprintf(fplog,"\n");		
    }
    fprintf(fplog,"%7c", ' '); 
    for(j = 0 ; j < n ; j ++) 
        fprintf(fplog, "%4d(%c) ", j, Ystr[j]);
    fprintf(fplog,"\n");
}/*}}}*/

int GetFastaSeqFromBuffer(char* buffer, vector <string> &idList, vector <string> &seqList, int &seq_type)/*{{{*/
{
    char *str = buffer; 
    int buffersize = strlen(buffer);
    char *pch_beg = NULL;
    char *pch_end = NULL;
    char *pch_seqbeg = NULL;
#ifdef DEBUG_READFASTA
    int cntSeq = 0;
    cout << "Reading sequences ..." << endl;
#endif
    while (1){
        pch_beg = strchr(str, '>'); 
        if (pch_beg == NULL){
            break;
        }
        pch_end = strstr(pch_beg+1, "\n>");
        if (pch_end == NULL){
            pch_end = buffer + strlen(buffer) - 1;
        }
        int sizerecord = pch_end - pch_beg + 1 ;  
        char *seqWithAnno = new char [sizerecord+1];
//        my_strcpy(seqWithAnno, pch_beg, sizerecord, sizerecord);
        my_strcpy(seqWithAnno, pch_beg, sizerecord);
        string seqid = "";
        GetIDFromAnnotationLine(seqWithAnno, seqid);

        pch_seqbeg = strchr(pch_beg+1,'\n');
        if (pch_seqbeg == NULL){
            cerr << "Fatal! fasta seqfile format error. Exit." << endl;
            exit(1);
        }
        int sizerawseq = pch_end - pch_seqbeg;
//        my_strcpy(seqWithAnno, pch_seqbeg+1, sizerawseq, sizerawseq);
        my_strcpy(seqWithAnno, pch_seqbeg+1, sizerawseq);

        char *seq = new char [sizerawseq+1];
        SpanExcluding(seqWithAnno,seq , WHITE_SPACE);
        idList.push_back(seqid);
        seqList.push_back(string(seq));
#ifdef DEBUG_READFASTA
        cntSeq += 1;
        //cout << ">" <<seqid<<endl;
        //cout << seq << endl;
        if (cntSeq % 10000 == 1){
            cout << cntSeq << "..." << flush;
            //cout << "position = " << str - buffer << endl;
        }
#endif
        delete [] seqWithAnno;
        delete [] seq;
        str = pch_end;
        if (str - buffer >= buffersize){
            break;
        }
    }
#ifdef DEBUG_READFASTA
        cout << "Finished." <<endl;
        cout << "Total sequences read in: " << cntSeq << endl;
#endif
    return 0;
}/*}}}*/
int ReadFasta(string &file,  vector <string> &idList, vector <string> &seqList, int &seq_type)/*{{{*/
{
    int filesize = GetFileSize(file.c_str());
#ifdef DEBUG_READFASTA
    cout << "file " << file << ". size = " << filesize << endl;
#endif 
    int MAX_FASTA_FILE_SIZE = 1*1024*1024*1024;
    if (filesize <= 0) {
        cerr << "Size of the file " << file 
            << " <= 0." 
            << endl;
        return -1;
    } else if (filesize > MAX_FASTA_FILE_SIZE){
        cerr << "Size of the file " << file 
            << " is over the limit (" 
            << MAX_FASTA_FILE_SIZE << ")."  << endl;
        return -1;
    }else{
        FILE *fpin = fopen(file.c_str(), "rb");
        if (fpin != NULL){
            char *buffer = new char [filesize +1];
            int nread = fread(buffer, sizeof (char), filesize, fpin);
#ifdef DEBUG_READFASTA
            cout << "buffer of size " << filesize+1 << " allocated." << endl;
            cout << "Read in " << nread << " bytes from file " << file << endl;
#endif 
            fclose(fpin);
            if (nread != filesize) {
                cerr << "Read file " << file <<
                    " failed. nread (" << nread << ") != filesize (" <<
                    filesize << ")." << endl;
                return -1;
            }
            buffer[filesize] = '\0';
            GetFastaSeqFromBuffer( buffer, idList, seqList, seq_type);
            delete [] buffer ; 
            return idList.size();
        }else{
            cerr << "Failed to read fasta file " << file << endl;
            return -1;
        }
    }
}/*}}}*/


void WriteAlignmentTableLine(string &title1, string &title2, int iAlign, int seqLength1, int seqLength2, AlignFactor &alignFactor, FILE*fpout)/*{{{*/
{
    fprintf(fpout,"%-16s %-15s %6.1f %6.1f %9d %6d %6d %9.1f %6d %6d %6d %6.1f %6.1f\n",
            title1.c_str(),title2.c_str(),
            alignFactor.identity*100.0,
            alignFactor.similarity*100.0,
            iAlign,
            seqLength1,
            seqLength2,
            alignFactor.score,
            alignFactor.idt_cnt,
            alignFactor.sim_cnt,
            alignFactor.gap_cnt,
            alignFactor.idt_cnt*100.0/min(seqLength1, seqLength2),
            alignFactor.idt_cnt*100.0/(iAlign - alignFactor.gap_cnt) );
}/*}}}*/
void WriteAlignmentHeaderNeedle(float gapOpen, float gapExt, int iAlign, AlignFactor *alignFactor, const char *title1, const char *title2, int seqLength1, int seqLength2, FILE *fpout)/*{{{*/
{
    fprintf(fpout,"#=======================================\n");
    fprintf(fpout,"#\n");
    fprintf(fpout,"# Aligned_sequences: %d\n", 2);
    fprintf(fpout,"# 1: %s\n", title1);
    fprintf(fpout,"# 2: %s\n", title2);
    fprintf(fpout,"# Matrix: %s\n", matrixFile.c_str());
    fprintf(fpout,"# Gap_penalty: %.1f\n", gapOpen);
    fprintf(fpout,"# Extend_penalty: %.1f\n", gapExt);
    fprintf(fpout,"#\n");
    fprintf(fpout,"# Length: %d\n", iAlign);
    fprintf(fpout,"# Identity:    %4d/%d (%4.1f%%)\n", alignFactor->idt_cnt, iAlign, alignFactor->identity*100.0);
    fprintf(fpout,"# Similarity:  %4d/%d (%4.1f%%)\n", alignFactor->sim_cnt, iAlign, alignFactor->similarity*100.0);
    fprintf(fpout,"# Gaps:        %4d/%d (%4.1f%%)\n", alignFactor->gap_cnt, iAlign, alignFactor->gapPercent*100.0);
    fprintf(fpout,"# Score: %.1f\n", alignFactor->score);
    fprintf(fpout,"# \n");
    fprintf(fpout,"#\n");
    fprintf(fpout,"#=======================================\n");
    fprintf(fpout,"\n");
}/*}}}*/
void WriteAlignmentNeedle(const char *title1, const char *title2, char *aXstr, char *aYstr, int *aRel, int length,int lineLength, FILE *fpout)/*{{{*/
{
    int i;
    int b1 = 0, e1  = 0; // begin and end count for sequence 1, without gaps
    int b2 = 0, e2  = 0; // begin and end count for sequence 2, without gaps
    int cntNonGapAA = 0;
    int lineCnt = 0 ;
    int lineBegin = 0 , lineEnd = 0 ;
    Array1D<char>linestr_1darray(lineLength+10);
    char* linestr= linestr_1darray.array1D;
    while(lineEnd < length)
    {	
        lineBegin = lineCnt * lineLength ;
        lineEnd = min((lineCnt + 1) * lineLength,length);
        //-------------------------------------------- line 1, alignseq1
        my_strcpy(linestr, aXstr+lineBegin,  lineLength);
        cntNonGapAA = 0;
        for(i = lineBegin ;i < lineEnd ; i ++ ) 
        {
            if(aXstr[i] != CHAR_GAP) { cntNonGapAA ++; }
        }
        if(cntNonGapAA>0) { b1 +=1; }
        e1 += cntNonGapAA;
        fprintf(fpout,"%-13.13s %6d %s %6d\n",title1,b1, linestr,e1); /*21 chars*/
        b1 = e1;
        //-------------------------------------------- line 2, relations
        fprintf(fpout,"%-13.13s %6s ","","");
        for(i = lineBegin ;i < lineEnd ; i ++ ) 
        {
            if(aRel[i] == IDT)
            {
                fprintf(fpout,"%c",'|');
            }
            else if(aRel[i] == SIM)
            {
                fprintf(fpout,"%c",':');
            }
            else if (aRel[i] != GAP)
            {
                fprintf(fpout,"%c",'.');
            }
            else
            {
                fprintf(fpout,"%c",' ');
            }
        }
        fprintf(fpout,"\n");
        //-------------------------------------------- line 3, alignseq2
        my_strcpy(linestr, aYstr+lineBegin,  lineLength);
        cntNonGapAA = 0;
        for(i = lineBegin ;i < lineEnd ; i ++ ) 
        {
            if(aYstr[i] != CHAR_GAP) { cntNonGapAA ++; }
        }
        if(cntNonGapAA>0) { b2 +=1; }
        e2 += cntNonGapAA;
        fprintf(fpout,"%-13.13s %6d %s %6d\n",title2,b2, linestr,e2); /*21 chars*/
        b2 = e2;
        //--------------------------------------------
        fprintf(fpout,"\n");
        lineCnt ++ ;
    }
    fprintf(fpout,"\n");
}/*}}}*/
void WriteAlignmentInFasta(string &title1, string &title2, char *alignstr, float gapopen, float gapextend, int iAlign, int seqLength, AlignFactor *pAlignFactor, FILE *fpout)/*{{{*/
{
    if (fpout != NULL){
        fprintf(fpout, ">%s aligned_to=%s seqIDT=%.1f seqIDT1=%.1f seqIDT2=%.1f seqSIM=%.1f alignLen=%d seqLen=%d go=%.2g ge=%.2g\n",
                title1.c_str(), title2.c_str(), 
                pAlignFactor->identity*100, 
                pAlignFactor->identity_short*100,
                pAlignFactor->idt_cnt*100.0/(iAlign - pAlignFactor->gap_cnt),
                pAlignFactor->similarity*100,
                iAlign, seqLength, gapopen, gapextend);
        fprintf(fpout, "%s\n", alignstr);
    }
}/*}}}*/


int embAlignWalkNWMatrixUsingCompass(const char* p, const char* q,/*{{{*/
                                 char *m, char *n,
                                 ajuint lena, ajuint lenb,
                                 ajint *start1, ajint *start2,
                                 ajint const *compass, FILE *fpLog)
{
    ajint xpos = *start2;
    ajint ypos = *start1;
    ajint i;
    ajint j;
    ajuint cursor;

    strcpy(m, "");
    strcpy(n, "");
    int cntM = 0;
    int cntN = 0;
    
    for (i=lenb-1; i>xpos;)
    {
        n[cntN]=q[i--] ;
        m[cntM]=CHAR_GAP;
        cntN ++;
        cntM ++;
    }

    for (j=lena-1; j>ypos; )
    {
        m[cntM]=p[j--] ;
        n[cntN]=CHAR_GAP;
        cntN ++;
        cntM ++;
    }

    while (xpos >= 0 && ypos >= 0)
    {
        cursor = ypos * lenb + xpos;
        if(!compass[cursor]) /* diagonal */
        {
            m[cntM] = p[ypos--];
            n[cntN] = q[xpos--];
            cntN ++;
            cntM ++;
            continue;
        }
        else if(compass[cursor] == LEFT) /* Left, gap(s) in vertical */
        {
            m[cntM]= CHAR_GAP;
            n[cntN]= q[xpos--];
            cntN ++;
            cntM ++;
            continue;
        }
        else if(compass[cursor] == DOWN) /* Down, gap(s) in horizontal */
        {
            m[cntM]= p[ypos--];
            n[cntN]= CHAR_GAP;
            cntN ++;
            cntM ++;
            continue;
        } 
        else
        {
            fprintf(stderr, "Walk Error in NW");
        }
    }

    for (;xpos>=0;xpos--)
    {
        n[cntN]= q[xpos];
        m[cntM]= CHAR_GAP;
        cntN ++;
        cntM ++;
    }

    for (;ypos>=0;ypos--)
    {
        m[cntM]= p[ypos];
        n[cntN]= CHAR_GAP;
        cntN ++;
        cntM ++;
    }

    *start2 = xpos+1;
    *start1 = ypos+1;

    m[cntM] = '\0';
    n[cntN] = '\0';

    StrReverse(m); /* written with append, need to reverse */
    StrReverse(n);

    if (fpLog != NULL) {  
        fprintf(fpLog,"embAlignWalkNWMatrixUsingCompass:\n");
        fprintf(fpLog, "first sequence extended with gaps  (m): %s\n", m);
        fprintf(fpLog, "second sequence extended with gaps (n): %s\n", n);
    }
    
    return cntM;
}/*}}}*/
static float embAlignGetScoreNWMatrix(/*{{{*/
	const float *ix, const float *iy, const float *m,
        ajint lena, ajint lenb,
        ajint *start1, ajint *start2,
        AjBool endweight)
{
    ajint i,j, cursor;
    float score = INT_MIN;
    *start1 = lena-1;
    *start2 = lenb-1;

    if(endweight)
    {
        /* when using end gap penalties the score of the optimal global
         * alignment is stored in the final cell of the path matrix */
        cursor = lena * lenb - 1;
        if(m[cursor]>ix[cursor]&&m[cursor]>iy[cursor])
            score = m[cursor];
        else if(ix[cursor]>iy[cursor])
            score = ix[cursor];
        else
            score = iy[cursor];
    }
    else {

        for (i = 0; i < lenb; ++i)
        {
            cursor = (lena - 1) * lenb + i;
            if(m[cursor]>score)
            {
                *start2 = i;
                score = m[cursor];
            }
            if(ix[cursor]>score)
            {
                score = ix[cursor];
                *start2 = i;
            }
            if(iy[cursor]>score)
            {
                score = iy[cursor];
                *start2 = i;
            }
        }

        for (j = 0; j < lena; ++j)
        {
            cursor = j * lenb + lenb - 1;
            if(m[cursor]>score)
            {
                *start1 = j;
                *start2 = lenb-1;
                score = m[cursor];
            }
            if(ix[cursor]>score)
            {
                score = ix[cursor];
                *start1 = j;
                *start2 = lenb-1;
            }
            if(iy[cursor]>score)
            {
                score = iy[cursor];
                *start1 = j;
                *start2 = lenb-1;
            }
        }
    }
    return score;
}/*}}}*/
float embAlignPathCalcWithEndGapPenalties(const char *a, const char *b,/*{{{*/
                       ajint lena, ajint lenb,
                       float gapopen, float gapextend,
                       float endgapopen, float endgapextend,
                       ajint *start1, ajint *start2,
                       float **sub, 
                       float *m, float *ix, float *iy,
                       ajint *compass, AjBool show,
                       AjBool endweight, FILE *fpLog)
{
    ajint xpos;
    ajint ypos;
    ajint bconvcode;

    float match;
    float ixp;
    float iyp;
    float mp;
    ajint cursor;
    ajint cursorp;
    
    float testog;
    float testeg;
    float score;
    
    
    if (!endweight)
    {
        endgapopen=0;
        endgapextend=0;
        endweight=ajTrue;
    }

    int i;
    sizeAlphabet=strlen(alphabet);

    Array1D <int8> X_1darray(lena+1);
    Array1D <int8> Y_1darray(lenb+1);
    int8 *X         = X_1darray.array1D;         /* digital X sequence         */
    int8 *Y         = Y_1darray.array1D;         /* digital Y sequence         */

    for(i = 0 ; i < lena ; i ++) X[i] = Charcase2Digit(a[i], alphabet, sizeAlphabet);
    for(i = 0 ; i < lenb ; i ++) Y[i] = Charcase2Digit(b[i], alphabet, sizeAlphabet);

    match = sub[X[0]][Y[0]];

    ix[0] = -endgapopen-gapopen;
    iy[0] = -endgapopen-gapopen;
    m[0] = match;

    cursor =0;

    /* First initialise the first column */
    for (ypos = 1; ypos < lena; ++ypos)
    {
	match = sub[X[ypos]][Y[0]];
	cursor = ypos * lenb;
	cursorp = (ypos-1) * lenb;

	testog = m[cursorp] - gapopen;
	testeg = iy[cursorp] - gapextend;

	if(testog >= testeg)
	    iy[cursor] = testog;
	else
	    iy[cursor] = testeg;

	m[cursor] = match - (endgapopen + (ypos - 1) * endgapextend);
	ix[cursor] = -endgapopen - ypos * endgapextend - gapopen;
    }

    ix[cursor] -= endgapopen;
    ix[cursor] += gapopen;

    cursor=0;

    /* Now initialise the first row */
    for (xpos = 1; xpos < lenb; ++xpos)
    {
	match = sub[X[0]][Y[xpos]];
	cursor = xpos;
	cursorp = xpos -1;

	testog = m[cursorp] - gapopen;
        testeg = ix[cursorp] - gapextend;

        if(testog >= testeg)
            ix[cursor] = testog;
        else
            ix[cursor] = testeg;

	m[cursor] = match - (endgapopen + (xpos - 1) * endgapextend);
	iy[cursor] = -endgapopen - xpos * endgapextend -gapopen;
    }

    iy[cursor] -= endgapopen;
    iy[cursor] += gapopen;

    xpos = 1;

    /* Now construct match, ix, and iy matrices */
    while (xpos != lenb)
    {
        ypos = 1;
        bconvcode = *++Y;

        /* coordinates of the cells being processed */
        cursorp = xpos-1;
        cursor = xpos++;

        while (ypos < lena)
        {
            /* get match for current xpos/ypos */
            match = sub[X[ypos++]][bconvcode];

            cursor += lenb;

            /* match matrix calculations */
            mp = m[cursorp];
            ixp = ix[cursorp];
            iyp = iy[cursorp];

            if(mp > ixp && mp > iyp)
                m[cursor] = mp+match;                     
            else if(ixp > iyp)
                m[cursor] = ixp+match;
            else
                m[cursor] = iyp+match;

            /* iy matrix calculations */
            if(xpos==lenb)
            {
        	testog = m[++cursorp] - endgapopen;
        	testeg = iy[cursorp] - endgapextend;
            }
            else
            {
        	testog = m[++cursorp];
        	
        	if (testog<ix[cursorp])
        	    testog = ix[cursorp];
        	
        	testog -= gapopen;
        	testeg = iy[cursorp] - gapextend;
            }
            
            if(testog > testeg)
                iy[cursor] = testog;
            else
        	iy[cursor] = testeg;

            cursorp += lenb;

            /* ix matrix calculations */
            if(ypos==lena)
            {
        	testog = m[--cursorp] - endgapopen;
        	testeg = ix[cursorp] - endgapextend;
            }
            else
            {
        	testog = m[--cursorp];
        	
        	if (testog<iy[cursorp])
        	    testog = iy[cursorp];
        	
        	testog -= gapopen;
        	testeg = ix[cursorp] - gapextend;
            }
            
            if(testog > testeg )
                ix[cursor] = testog;
            else
        	ix[cursor] = testeg;

        }
    }

    score = embAlignGetScoreNWMatrix(ix, iy, m, lena, lenb, start1, start2, endweight);

    xpos = *start2;
    ypos = *start1;

    /* In the following loop the three matrices (m, ix, iy) are traced back
     * and path/alignment decision/selection is made.
     * 0 means match: go up and left in the matrix
     * 1 means: go left in the matrix, i.e. gap in the first sequence(seq a)
     * 2 means: go up in the matrix, i.e. gap in the second sequence(seq b)
     */
    cursorp=0;
    cursor=1;
    
    while (xpos>=0 && ypos>=0)
    {
	cursor = ypos*lenb+xpos;
	mp = m[cursor];

	if(cursorp == LEFT && E_FPEQ((ypos==0||(ypos==lena-1)?
		endgapextend:gapextend), (ix[cursor]-ix[cursor+1]),U_FEPS))
	{
	    compass[cursor] = LEFT;
	    xpos--;
	}
	else if(cursorp== DOWN && E_FPEQ((xpos==0||(xpos==lenb-1)?
		endgapextend:gapextend), (iy[cursor]-iy[cursor+lenb]),U_FEPS))
	{
	    compass[cursor] = DOWN;
	    ypos--;
	}
	else if(mp >= ix[cursor] && mp>= iy[cursor])
	{

	    if(cursorp == LEFT && E_FPEQ(mp,ix[cursor],U_FEPS))
	    {
		compass[cursor] = LEFT;
		xpos--;
	    }
	    else if(cursorp == DOWN && E_FPEQ(mp,iy[cursor],U_FEPS))
	    {
		compass[cursor] = DOWN;
		ypos--;
	    }
	    else
	    {
		compass[cursor] = 0;
		ypos--;
		xpos--;
	    }

	}
	else if(ix[cursor]>=iy[cursor] && xpos>-1)
	{
	    compass[cursor] = LEFT;
	    xpos--;
	}
	else if(ypos>-1)
	{
	    compass[cursor] = DOWN;
	    ypos--;
	}
	else
	{
	    fprintf(stderr,"something is seriously wrong in the traceback algorithm\n");
	    exit(1);
	}
	cursorp = compass[cursor];
    }

    if(show)
    {
        printPathMatrix(m, compass, a, b, lena, lenb, fpLog);
        printPathMatrix(ix, compass, a, b, lena, lenb, fpLog);
        printPathMatrix(iy, compass, a, b, lena, lenb, fpLog);
    }
    return score;
}/*}}}*/
float embAlignPathCalcWithEndGapPenalties2(const char *a, const char *b,/*{{{*/
                       int lena, int lenb,
                       float* gapopenArray1, float *gapopenArray2, float gapextend,
                       float endgapopen, float endgapextend,
                       int *start1, int *start2,
                       float * const *sub, 
                       float *m, float *ix, float *iy,
                       int *compass, int show,
                       AjBool endweight, FILE *fpLog)
{
    ajint xpos;
    ajint ypos;
    ajint bconvcode;

    float match;
    float ixp;
    float iyp;
    float mp;
    ajint cursor;
    ajint cursorp;

    float testog;
    float testeg;
    float score;


    if (!endweight)
    {
        endgapopen=0;
        endgapextend=0;
        endweight=true;
    }

    int i;
    sizeAlphabet=strlen(alphabet);

    Array1D <int8> X_1darray(lena+1);
    Array1D <int8> Y_1darray(lenb+1);
    int8 *X         = X_1darray.array1D;         /* digital X sequence         */
    int8 *Y         = Y_1darray.array1D;         /* digital Y sequence         */

    for(i = 0 ; i < lena ; i ++) X[i] = Charcase2Digit(a[i], alphabet, sizeAlphabet);
    for(i = 0 ; i < lenb ; i ++) Y[i] = Charcase2Digit(b[i], alphabet, sizeAlphabet);

    match = sub[X[0]][Y[0]];

    /*ix[0] = -endgapopen-gapopen;*/
    /*iy[0] = -endgapopen-gapopen;*/
    ix[0] = -endgapopen-gapopenArray1[0];
    iy[0] = -endgapopen-gapopenArray2[0];
    m[0] = match;

    cursor =0;

    /* First initialise the first column */
    for (ypos = 1; ypos < lena; ++ypos)
    {
        match = sub[X[ypos]][Y[0]];
        cursor = ypos * lenb;
        cursorp = (ypos-1) * lenb;

        /*testog = m[cursorp] - gapopen;*/
        testog = m[cursorp] - gapopenArray1[ypos];
        testeg = iy[cursorp] - gapextend;

        if(testog >= testeg)
            iy[cursor] = testog;
        else
            iy[cursor] = testeg;

        m[cursor] = match - (endgapopen + (ypos - 1) * endgapextend);
        /*ix[cursor] = -endgapopen - ypos * endgapextend - gapopen;*/
        ix[cursor] = -endgapopen - ypos * endgapextend - gapopenArray1[ypos];
    }

    ix[cursor] -= endgapopen;
    /*ix[cursor] += gapopen;*/
    ix[cursor] += gapopenArray1[ypos-1];

    cursor=0;

    /* Now initialise the first row */
    for (xpos = 1; xpos < lenb; ++xpos)
    {
        match = sub[X[0]][Y[xpos]];
        cursor = xpos;
        cursorp = xpos -1;

        /*testog = m[cursorp] - gapopen;*/
        testog = m[cursorp] - gapopenArray2[xpos];
        testeg = ix[cursorp] - gapextend;

        if(testog >= testeg)
            ix[cursor] = testog;
        else
            ix[cursor] = testeg;

        m[cursor] = match - (endgapopen + (xpos - 1) * endgapextend);
        /*iy[cursor] = -endgapopen - xpos * endgapextend -gapopen;*/
        iy[cursor] = -endgapopen - xpos * endgapextend -gapopenArray2[xpos];
    }

    iy[cursor] -= endgapopen;
    /*iy[cursor] += gapopen;*/
    iy[cursor] += gapopenArray2[xpos-1];

    xpos = 1;

    /* Now construct match, ix, and iy matrices */
    while (xpos != lenb)
    {
        ypos = 1;
        bconvcode = *++Y;

        /* coordinates of the cells being processed */
        cursorp = xpos-1;
        cursor = xpos++;

        while (ypos < lena)
        {
            /* get match for current xpos/ypos */
            match = sub[X[ypos++]][bconvcode];

            cursor += lenb;

            /* match matrix calculations */
            mp = m[cursorp];
            ixp = ix[cursorp];
            iyp = iy[cursorp];

            if(mp > ixp && mp > iyp)
                m[cursor] = mp+match;                     
            else if(ixp > iyp)
                m[cursor] = ixp+match;
            else
                m[cursor] = iyp+match;

            /* iy matrix calculations */
            if(xpos==lenb)
            {
                testog = m[++cursorp] - endgapopen;
                testeg = iy[cursorp] - endgapextend;
            }
            else
            {
                testog = m[++cursorp];

                if (testog<ix[cursorp])
                    testog = ix[cursorp];

                /*testog -= gapopen;*/
                testog -= gapopenArray2[xpos];
                testeg = iy[cursorp] - gapextend;
            }

            if(testog > testeg)
                iy[cursor] = testog;
            else
                iy[cursor] = testeg;

            cursorp += lenb;

            /* ix matrix calculations */
            if(ypos==lena)
            {
                testog = m[--cursorp] - endgapopen;
                testeg = ix[cursorp] - endgapextend;
            }
            else
            {
                testog = m[--cursorp];

                if (testog<iy[cursorp])
                    testog = iy[cursorp];

                /*testog -= gapopen;*/
                testog -= gapopenArray1[ypos];
                testeg = ix[cursorp] - gapextend;
            }

            if(testog > testeg )
                ix[cursor] = testog;
            else
                ix[cursor] = testeg;

        }
    }

    score = embAlignGetScoreNWMatrix(ix, iy, m, lena, lenb, start1, start2, endweight);

    xpos = *start2;
    ypos = *start1;

    /* In the following loop the three matrices (m, ix, iy) are traced back
     * and path/alignment decision/selection is made.
     * 0 means match: go up and left in the matrix
     * 1 means: go left in the matrix, i.e. gap in the first sequence(seq a)
     * 2 means: go up in the matrix, i.e. gap in the second sequence(seq b)
     */
    cursorp=0;
    cursor=1;

    while (xpos>=0 && ypos>=0)
    {
        cursor = ypos*lenb+xpos;
        mp = m[cursor];

        if(cursorp == LEFT && E_FPEQ((ypos==0||(ypos==lena-1)?  endgapextend:gapextend), (ix[cursor]-ix[cursor+1]),U_FEPS))
        {
            compass[cursor] = LEFT;
            xpos--;
        }
        else if(cursorp== DOWN && E_FPEQ((xpos==0||(xpos==lenb-1)?  endgapextend:gapextend), (iy[cursor]-iy[cursor+lenb]),U_FEPS))
        {
            compass[cursor] = DOWN;
            ypos--;
        }
        else if(mp >= ix[cursor] && mp>= iy[cursor])
        {

            if(cursorp == LEFT && E_FPEQ(mp,ix[cursor],U_FEPS))
            {
                compass[cursor] = LEFT;
                xpos--;
            }
            else if(cursorp == DOWN && E_FPEQ(mp,iy[cursor],U_FEPS))
            {
                compass[cursor] = DOWN;
                ypos--;
            }
            else
            {
                compass[cursor] = 0;
                ypos--;
                xpos--;
            }

        }
        else if(ix[cursor]>=iy[cursor] && xpos>-1)
        {
            compass[cursor] = LEFT;
            xpos--;
        }
        else if(ypos>-1)
        {
            compass[cursor] = DOWN;
            ypos--;
        }
        else
        {
            fprintf(stderr, "something is seriously wrong in the traceback algorithm\n");
            exit(1);
        }
        cursorp = compass[cursor];
    }

    if(show)
    {
        printPathMatrix(m, compass, a, b, lena, lenb, fpLog);
        printPathMatrix(ix, compass, a, b, lena, lenb, fpLog);
        printPathMatrix(iy, compass, a, b, lena, lenb, fpLog);
    }
    return score;
}/*}}}*/

int Alignment_MyNeedle(const char *p, const char *q, const char *alphabet, int lena, int lenb, const char *title1, const char *title2, char *alnP,char *alnQ,int *alignRel, AlignFactor &alignFactor, float gapopen, float gapextend, float endgapopen, float endgapextend, float **sub, int seq_type, FILE *fpLog)/*{{{*/
/*constant gap penalties*/
{
    float score = 0.0;
    int start1;
    int start2;

    int maxSize= (lena+1)*(lenb+1)+1;
    Array1D <int> compass_1darray(maxSize);
    Array1D <float> ix_1darray(maxSize);
    Array1D <float> iy_1darray(maxSize);
    Array1D <float> m_1darray(maxSize);

    int *compass = compass_1darray.array1D;
    float* ix = ix_1darray.array1D;
    float* iy = iy_1darray.array1D;
    float* m = m_1darray.array1D;

    score = embAlignPathCalcWithEndGapPenalties(p, q, lena, lenb, gapopen,  gapextend, endgapopen, endgapextend, &start1, &start2, sub, m, ix, iy, compass, show, endweight, fpLog);

    int alignLength = embAlignWalkNWMatrixUsingCompass(p, q, alnP, alnQ, lena, lenb, &start1, &start2, compass, fpLog);
    AlignFactor *pAlignFactor = &alignFactor;

    AlignAna(alnP,alnQ,alphabet, alignLength ,sub,pAlignFactor, alignRel);
    pAlignFactor->identity_short = pAlignFactor->identity * alignLength / min(lena,lenb) ;
    pAlignFactor->similarity_short = pAlignFactor->similarity * alignLength / min(lena,lenb) ;
    pAlignFactor->score = score;
    return alignLength ;
}/*}}}*/
int Alignment_MyNeedle2(const char *p, const char *q, const char *alphabet, int lena, int lenb, const char *title1, const char *title2, char *alnP,char *alnQ,int *alignRel, AlignFactor &alignFactor, float *gapopenArray1, float *gapopenArray2, float gapextend, float endgapopen, float endgapextend, float **sub, int seq_type  , FILE *fpLog)/*{{{*/
/*position specific gap penalties*/
{
    float score = 0.0;
    int start1;
    int start2;

    int maxSize= (lena+1)*(lenb+1)+1;
    Array1D <int> compass_1darray(maxSize);
    Array1D <float> ix_1darray(maxSize);
    Array1D <float> iy_1darray(maxSize);
    Array1D <float> m_1darray(maxSize);

    int *compass = compass_1darray.array1D;
    float* ix = ix_1darray.array1D;
    float* iy = iy_1darray.array1D;
    float* m = m_1darray.array1D;

    score = embAlignPathCalcWithEndGapPenalties2(p, q, lena, lenb, gapopenArray1,gapopenArray2,  gapextend, endgapopen, endgapextend, &start1, &start2, sub, m, ix, iy, compass, show, endweight, fpLog);

    int alignLength = embAlignWalkNWMatrixUsingCompass(p, q, alnP, alnQ, lena, lenb, &start1, &start2, compass, fpLog);
    AlignFactor *pAlignFactor = &alignFactor;

    AlignAna(alnP,alnQ,alphabet, alignLength ,sub,pAlignFactor, alignRel);
    pAlignFactor->identity_short = pAlignFactor->identity * alignLength / min(lena,lenb) ;
    pAlignFactor->similarity_short = pAlignFactor->similarity * alignLength / min(lena,lenb) ;
    pAlignFactor->score = score;
    return alignLength ;
}/*}}}*/
int GetSeqFromBuffer(char*buffer, string &aaSeq)/*{{{*/
{
#ifdef DEBUG_GETSEQFROMBUFFER
    cout << "buffer=\"" << buffer << "\"" <<endl;
#endif
    char *pch = buffer;
    if (buffer[0] == '>'){
        pch = strchr(buffer, '\n');
        pch ++;
    }
    int lengthrawseq = strlen(pch);
    char *tmpstr= new char [lengthrawseq+1];
    int cnt= 0;
    while (*pch){
        if (*pch >= 'A' && *pch <= 'z'){
            tmpstr[cnt] = *pch;
            cnt ++;
        }
        pch ++;
    }
    tmpstr[cnt] = '\0';
    aaSeq = tmpstr;
    delete [] tmpstr;
    return 0;
}/*}}}*/
int GetAASeqFromDatabase(string &id, string &aaSeq, map<string, dbindex>&dbindexmap, vector<FILE*> fpList)/*{{{*/
{
    if (dbindexmap.find(id) == dbindexmap.end()){
        fprintf(stderr, "id %s can not find in the seq database. Ignore.\n", id.c_str());
        return -1;
    }
    int status_fseek;
    dbindex *pDBindex = &dbindexmap[id];
#ifdef DEBUG_DBINDEX
    cout << id << " dbfileindex = " << pDBindex->dbfileindex << " offset=" << pDBindex->offset << endl; 
#endif
    if ((status_fseek = fseek(fpList[pDBindex->dbfileindex], pDBindex->offset, SEEK_SET)) == 0) {
        size_t nread;
        char *buffer = new char [ pDBindex->size+1];
        nread = fread(buffer, sizeof(char), pDBindex->size,
                fpList[pDBindex->dbfileindex]);
        buffer[pDBindex->size] = '\0';
        if (nread != pDBindex->size){
            cerr <<"fread failed" << endl;
            exit(1);
        }
        GetSeqFromBuffer(buffer, aaSeq);
        delete [] buffer;
    }else {
        cerr <<"Fatal! fseek failed for id \'"
            << id
            << "\' in the database. Exit."
            << endl;
        exit(1);
    }
    return 0;
}/*}}}*/
int GetAASeqFromDatabase2(string &aaSeq, map<string, dbindex>:: iterator itDBIndexMap, vector<FILE*> fpList)/*{{{*/
{
    int status_fseek;
    dbindex *pDBindex = &(itDBIndexMap->second);
#ifdef DEBUG_DBINDEX
    cout << id << " dbfileindex = " << pDBindex->dbfileindex << " offset=" << pDBindex->offset << endl; 
#endif
    if ((status_fseek = fseek(fpList[pDBindex->dbfileindex], pDBindex->offset, SEEK_SET)) == 0) {
        size_t nread;
        char *buffer = new char [ pDBindex->size+1];
        nread = fread(buffer, sizeof(char), pDBindex->size,
                fpList[pDBindex->dbfileindex]);
        buffer[pDBindex->size] = '\0';
        if (nread != pDBindex->size){
            cerr <<"fread failed" << endl;
            exit(1);
        }
        GetSeqFromBuffer(buffer, aaSeq);
        delete [] buffer;
    }else {
        cerr <<"Fatal! fseek failed for id \'"
            << itDBIndexMap->first
            << "\' in the database. Exit."
            << endl;
        exit(1);
    }
    return 0;
}/*}}}*/

int ReadDatabaseIndex_text(string &indexfile, map<string, dbindex> &dbindexmap,int &maxdbfileindex )/*{{{*/
{
    /*read in the whole file first*/
    ifstream ifp (indexfile.c_str(), ios::binary);
    dbindex tmpdbindex;
    if (ifp.is_open()){
        // get length of file:
        ifp.seekg (0, ios::end);
        int length = ifp.tellg();
        ifp.seekg (0, ios::beg);
        char *buffer = new char [ length+1];
        ifp.read(buffer,length);
        ifp.close();
        buffer[length] = '\0';
        char *pch;
        pch = strtok(buffer, "\n");
        while (pch != NULL){
            if ( strlen(pch)> 0 && strncmp(pch, "DEF", 3) != 0){
                char *id = new char [ strlen(pch)] ; 
                if (sscanf( pch, "%s %d %ld %ld", id, &(tmpdbindex.dbfileindex)
                        ,&(tmpdbindex.offset),&(tmpdbindex.size)) != 4){
                    cerr << "Failed to read the indexfile "
                        << indexfile << endl;
                    return -1;
                } else {
                    dbindexmap.insert(pair<string, dbindex>(string(id), tmpdbindex));
                    if (tmpdbindex.dbfileindex > maxdbfileindex){
                        maxdbfileindex = tmpdbindex.dbfileindex;
                    }
                }
                delete [] id;
            }
            pch = strtok(NULL, "\n");
        }
        delete [] buffer;
        return dbindexmap.size();
    } else {
        cerr << "Failed to open database index file "
            << indexfile 
            << endl;
        return -1;
    }
}/*}}}*/
int ReadDatabaseIndex_binary(string &indexfile, map<string, dbindex> &dbindexmap,int &maxdbfileindex )/*{{{*/
{
    /*read in the whole file first . binary index*/
    ifstream ifp (indexfile.c_str(), ios::binary);
    if (ifp.is_open()){
        // get length of file:
        ifp.seekg (0, ios::end);
        int length = ifp.tellg();
        ifp.seekg (0, ios::beg);
        BYTE *buffer = new BYTE [ length+1];
        ifp.read(reinterpret_cast<char*>(buffer),length);
        buffer[length] = '\0';
        ifp.close();

        BYTE *pBuffer = buffer;
        unsigned int sizedumpedtext=0;
        memcpy(&sizedumpedtext,pBuffer,sizeof(unsigned int));
        //memcpy(dumpedtext,pBuffer,sizeof(unsigned int));
        pBuffer+=sizeof(unsigned int) + sizedumpedtext*sizeof(char);

        unsigned int dumpedidlistsize=0;
        memcpy(&dumpedidlistsize,pBuffer,sizeof(unsigned int));
        pBuffer+=sizeof(unsigned int);
        char *dumpedidlist  = new char[dumpedidlistsize+1];
        memcpy(dumpedidlist,pBuffer,dumpedidlistsize);
        pBuffer+=sizeof(char)*dumpedidlistsize;

        vector <string> idList;
        char *pch = NULL;
        pch = strtok(dumpedidlist, "\n");
        while (pch != NULL){
            if (strlen(pch) > 0){
                idList.push_back(pch);
            }
            pch = strtok(NULL, "\n");
        }
        delete [] dumpedidlist;

        unsigned int numRecord=0;
        memcpy(&numRecord,pBuffer,sizeof(unsigned int));
        pBuffer+=sizeof(unsigned int);
        if (numRecord != idList.size()){
            cerr << "numRecord("
                << numRecord
                << ") != idList.size ("
                << idList.size()
                << ")."
                << endl;
        }

        unsigned char *arrayFileIndex=new unsigned char[numRecord];
        unsigned int  *arrayoffset =new unsigned int[numRecord];
        unsigned int  *arraysize =new unsigned int[numRecord];
        memcpy(arrayFileIndex,pBuffer,sizeof(unsigned char)*numRecord);
        pBuffer+=sizeof(unsigned char)*numRecord;
        memcpy(arrayoffset,pBuffer,sizeof(unsigned int)*numRecord);
        pBuffer+=sizeof(unsigned int)*numRecord;
        memcpy(arraysize,pBuffer,sizeof(unsigned int)*numRecord);
        pBuffer+=sizeof(unsigned int)*numRecord;

        if ((pBuffer-buffer) != length ){
            cerr << "Failed to read binary index file "
                << " filesize = " << length
                << " readsize = " << pBuffer-buffer
                << endl;
        }
        dbindex tmpdbindex;
        for (unsigned int i = 0; i < numRecord; i++) {
            tmpdbindex.dbfileindex=arrayFileIndex[i];
            tmpdbindex.offset=arrayoffset[i];
            tmpdbindex.size=arraysize[i];
            dbindexmap.insert(pair<string, dbindex>(idList[i], tmpdbindex));
            if (tmpdbindex.dbfileindex > maxdbfileindex){
                maxdbfileindex = tmpdbindex.dbfileindex;
            }
        }
        delete [] arrayFileIndex;
        delete [] arraysize;
        delete [] arrayoffset;
        delete [] buffer;
        return dbindexmap.size();
    } else {
        cerr << "Failed to open database index file "
            << indexfile 
            << endl;
        return -1;
    }
}/*}}}*/
int ReadDatabaseIndex(string &dbname, map<string,dbindex>&dbindexmap, int &maxDBIndexNumber)/*{{{*/
{
    string indexfile_binary=dbname+string(".indexbin");
    string indexfile_text=dbname+string(".index");
    if ( IsFileExist( indexfile_binary.c_str()) ) {
        return  ReadDatabaseIndex_binary(indexfile_binary, dbindexmap, maxDBIndexNumber); 
    }else if ( IsFileExist( indexfile_text.c_str()) ) {
        return  ReadDatabaseIndex_text(indexfile_text, dbindexmap, maxDBIndexNumber); 
    }else{
        cerr << "Index file for seqDatabase "
            << dbname
            << "does not exist." << endl;
        return -1;
    }

}/*}}}*/

int GetDBFPList( vector <FILE*> &fpList, string dbname, int maxdbfileindex)/*{{{*/
{
    for (int i=0;i<=maxdbfileindex;i++){
        string filename=dbname + int2string(i) + string(".db");
        FILE *fp = fopen(filename.c_str(),"rb");
        if (fp != NULL){
            fpList.push_back(fp);
        }else{
            cerr << "Failed to open dbfile "
                << filename
                << ". Exit."
                << endl;
            exit(1);
        }
    }
    return 0;
}/*}}}*/
int GetIDSeqMap(map <string, dbindex> &dbindexmap,vector <FILE*>&fpList, map <string,string> &idseqmap)/*{{{*/
{
    map <string, dbindex> :: iterator it;
    string aaSeq;
    for (it = dbindexmap.begin(); it != dbindexmap.end(); it++) {
        if (GetAASeqFromDatabase2(aaSeq, it, fpList) == 0 ){
            idseqmap.insert(pair<string, string>(it->first, aaSeq));
        }
    }
    return idseqmap.size();
}/*}}}*/
int AlignAllPair (vector<Pair> &pairList, map <string, dbindex> &dbindexmap,vector<FILE*>& fpList, map <string, string>&idseqmap, int dbfilesize, float ** SM, FILE *fpout, FILE *fpLog , FILE *fpTable)/*{{{*/
{
    int i;
    int numPair = pairList.size();
    int sizeTable = int(pow(double(23), double(wordsize)));
    int *freqTable1 = new int [sizeTable+1];
    int *freqTable2 = new int [sizeTable+1];

    for (i=0;i <numPair; i++){
        string id1 = pairList[i].mem1;
        string id2 = pairList[i].mem2;
#ifdef DEBUG
        cout << "align " << id1 << " " << id2 << endl;
#endif
        string aaSeq1;
        string aaSeq2;
        if (dbfilesize < THRESHOLD_DB_FILESIZE){
            if ((aaSeq1 = idseqmap[ id1 ]) == "" ) {
                continue;
            }
            if ((aaSeq2 = idseqmap[ id2 ]) == "" ) {
                continue;
            }
        } else {
            if (GetAASeqFromDatabase(id1, aaSeq1, dbindexmap, fpList) == -1 ){
                continue;
            }
            if (GetAASeqFromDatabase(id2, aaSeq2, dbindexmap, fpList) == -1 ){
                continue;
            }
        }
        string title1 = id1;
        string title2 = id2;
        bool isHasGapOpenArray1 = false;
        bool isHasGapOpenArray2 = false;
        float *gapopenArray1 = NULL;
        float *gapopenArray2 = NULL;
        int seqLength1 = aaSeq1.size();
        int seqLength2 = aaSeq2.size();
        if(seqLength1 * seqLength2 > MAX_ALIGN_SIZE) {
            fprintf(stderr,"%s (%d aa) - %s (%d aa), length over limit. Ignore.\n", id1.c_str(), seqLength1, id2.c_str(), seqLength2);
            return (-1);
        }
        float kmerBitScore  = 0.0;
        if (isPrintKMerBitScore){
            //kmerBitScore = KMerVectorPairwiseComparison(aaSeq1, aaSeq2, id1, id2, wordsize, fpLog);
            //kmerBitScore = KMerVectorPairwiseComparison_c(aaSeq1, aaSeq2, id1, id2, wordsize, fpLog);
            kmerBitScore = KMerVectorPairwiseComparison_table(aaSeq1, aaSeq2, id1, id2, freqTable1, freqTable2, sizeTable, wordsize, fpLog);
        }

        Array1D <char> alignXstr_1darray(seqLength1+seqLength2+1);
        Array1D <char> alignYstr_1darray(seqLength1+seqLength2+1);
        Array1D <int>  alignRel_1darray(seqLength1+seqLength2+1);
        char *alignXstr = alignXstr_1darray.array1D;
        char *alignYstr = alignYstr_1darray.array1D;
        int  *alignRel  = alignRel_1darray.array1D;
        AlignFactor alignFactor;
        InitAlignFactor(&alignFactor);

#ifdef DEBUG
        cout << "before alignment" << endl;
#endif
        int iAlign = 0;/*alignment length*/
        if ((noGapOpenArray == false)&&(isHasGapOpenArray1 || isHasGapOpenArray2)) {
            iAlign = Alignment_MyNeedle2(aaSeq1.c_str(), aaSeq2.c_str(), alphabet, seqLength1, seqLength2, title1.c_str(), title2.c_str(), alignXstr,alignYstr,alignRel, alignFactor, gapopenArray1, gapopenArray2, gapextend, endgapopen, endgapextend, SM, seq_type, fpLog);
        } else {
#ifdef DEBUG
            cout << "Run Alignment_MyNeedle" << endl;
#endif
            iAlign = Alignment_MyNeedle(aaSeq1.c_str(), aaSeq2.c_str(), alphabet, seqLength1, seqLength2, title1.c_str(), title2.c_str(), alignXstr,alignYstr,alignRel, alignFactor, gapopen, gapextend, endgapopen, endgapextend, SM, seq_type, fpLog);
        }
#ifdef DEBUG
        cout << " Alignment finished" << endl;
#endif
        if (isPrintKMerBitScore){
            fprintf(stdout,"%-12s - %12s vector bitScore = %10.3g seqidt = %10.3g len1 = %d len2 = %d iAlign = %d\n", id1.c_str(),id2.c_str(), kmerBitScore, alignFactor.identity*100, seqLength1, seqLength2, iAlign);
        }
        if (outformat == 0){
                WriteAlignmentHeaderNeedle(gapopen, gapextend, iAlign, &alignFactor, title1.c_str(), title2.c_str(), seqLength1, seqLength2, fpout);
                WriteAlignmentNeedle(title1.c_str(), title2.c_str(), alignXstr, alignYstr, alignRel, iAlign,nchar, fpout);
        }else if (outformat==1){
            WriteAlignmentInFasta( title1, title2, alignXstr, gapopen, gapextend, iAlign, seqLength1, &alignFactor, fpout  );
            WriteAlignmentInFasta( title2, title1, alignYstr, gapopen, gapextend, iAlign, seqLength2, &alignFactor, fpout  );
        }
        if (fpTable != NULL){
            WriteAlignmentTableLine(title1, title2, iAlign, seqLength1, seqLength2, alignFactor, fpTable);
        }
    }
    delete [] freqTable1;
    delete [] freqTable2;
    return 0;
}/*}}}*/
int AlignTwoSeqFile(string &file1, string &file2, string &title1, string &title2, float ** SM, FILE *fpout, FILE* fpLog, FILE* fpTable)/*{{{*/
{
    int seq_type1;
    int seq_type2;
    int seqLength1;
    int seqLength2;
    Array1D <char> aaSeq1_1darray(LONGEST_SEQ+1);
    Array1D <char> aaSeq2_1darray(LONGEST_SEQ+1);
    char *aaSeq1 = aaSeq1_1darray.array1D;
    char *aaSeq2 = aaSeq2_1darray.array1D;
    int maxSizeAnnoLine = 500;
    Array1D <char> annoLine1_1darray(maxSizeAnnoLine+1);
    Array1D <char> annoLine2_1darray(maxSizeAnnoLine+1);
    char *annoLine1 = annoLine1_1darray.array1D;
    char *annoLine2 = annoLine2_1darray.array1D;

    Array1D <float> gapopenArray1_1darray(LONGEST_SEQ);
    Array1D <float> gpeArray1_1darray(LONGEST_SEQ);
    Array1D <float> tgpeArray1_1darray(LONGEST_SEQ);
    float *gapopenArray1 = gapopenArray1_1darray.array1D;
    float *gpeArray1 = gpeArray1_1darray.array1D;
    float *tgpeArray1 = tgpeArray1_1darray.array1D;
    bool isHasGapOpenArray1 = false;

    Array1D <float> gapopenArray2_1darray(LONGEST_SEQ);
    Array1D <float> gpeArray2_1darray(LONGEST_SEQ);
    Array1D <float> tgpeArray2_1darray(LONGEST_SEQ);
    float *gapopenArray2 = gapopenArray2_1darray.array1D;
    float *gpeArray2 = gpeArray2_1darray.array1D;
    float *tgpeArray2 = tgpeArray2_1darray.array1D;
    bool isHasGapOpenArray2 = false;

    Array1D <char> rtname1_1darray(file1.size()+1);
    Array1D <char> rtname2_1darray(file2.size()+1);
    char *rtname1 = rtname1_1darray.array1D;
    char *rtname2 = rtname2_1darray.array1D;
    rootname(file1.c_str(),rtname1,file1.size());
    rootname(file2.c_str(),rtname2,file2.size());

    seqLength1 = ReadSeq_FASTA(file1.c_str(), aaSeq1, gapopenArray1, gpeArray1, tgpeArray1, isHasGapOpenArray1, &seq_type1, LONGEST_SEQ, annoLine1, maxSizeAnnoLine);

    if (!isHasGapOpenArray1) {
        for (int jj = 0; jj < seqLength1; jj++) {
            gapopenArray1[jj] = gapopen;
        }
    }

    GetIDFromAnnotationLine(annoLine1, title1); 
    if (title1 ==  "" ) {
        title1 = rtname1;
    }

    if (fpLog != NULL && isHasGapOpenArray1) {
        fprintf(fpLog, "seq: %s has gapopen\n", title1.c_str());
        for (int kk = 0; kk < seqLength1; kk ++) {
            fprintf(fpLog, "gapopenArray1[%d] = %.1f\n",kk, gapopenArray1[kk] );
        }
    }

    int cntSeq2 = 0;
    FILE *fpseq2 = fopen(file2.c_str(), "r");
    if (fpseq2 == NULL){
        fprintf(stderr,"Failed to read seqfile2 %s. Ignore.", file2.c_str());
    }else{
        while ((seqLength2 = ReadNextSeq_FASTA(fpseq2, aaSeq2, gapopenArray2, gpeArray2, tgpeArray2, isHasGapOpenArray2, &seq_type2, LONGEST_SEQ, annoLine2, maxSizeAnnoLine))!= EOF) {
            cntSeq2++;
            if(seqLength1 * seqLength2 > MAX_ALIGN_SIZE) {
                fprintf(stderr,"%s (%d aa) - %s (%d aa), length over limit. Ignore.\n", file1.c_str(), seqLength1, file2.c_str(), seqLength2);
                continue;
            }
            GetIDFromAnnotationLine(annoLine2, title2); 
            if(title2 == ""){
                title2 = rtname2 +  string("_") +  int2string(cntSeq2);
            }

            if (!isHasGapOpenArray2) {
                for (int jj = 0; jj < seqLength2; jj++) {
                    gapopenArray2[jj] = gapopen;
                }
            }

            if (fpLog != NULL && isHasGapOpenArray2) {
                fprintf(fpLog, "seq: %s has gapopen\n", title2.c_str());
                for (int kk = 0; kk < seqLength2; kk ++) {
                    fprintf(fpLog, "gapopenArray2[%d] = %.1f\n",kk, gapopenArray2[kk] );
                }
            }

            Array1D <char> alignXstr_1darray(seqLength1+seqLength2+1);
            Array1D <char> alignYstr_1darray(seqLength1+seqLength2+1);
            Array1D <int>  alignRel_1darray(seqLength1+seqLength2+1);
            char *alignXstr = alignXstr_1darray.array1D;
            char *alignYstr = alignYstr_1darray.array1D;
            int  *alignRel  = alignRel_1darray.array1D;
            AlignFactor alignFactor;
            InitAlignFactor(&alignFactor);

            int iAlign = 0;/*alignment length*/
            if ((noGapOpenArray == false)&&(isHasGapOpenArray1 || isHasGapOpenArray2)) {
                iAlign = Alignment_MyNeedle2(aaSeq1, aaSeq2, alphabet, seqLength1, seqLength2, title1.c_str(), title2.c_str(), alignXstr,alignYstr,alignRel, alignFactor, gapopenArray1, gapopenArray2, gapextend, endgapopen, endgapextend, SM, seq_type, fpLog);
            } else {
                iAlign = Alignment_MyNeedle(aaSeq1, aaSeq2, alphabet, seqLength1, seqLength2, title1.c_str(), title2.c_str(), alignXstr,alignYstr,alignRel, alignFactor, gapopen, gapextend, endgapopen, endgapextend, SM, seq_type, fpLog);
            }
            if (outformat == 0){
                WriteAlignmentHeaderNeedle(gapopen, gapextend, iAlign, &alignFactor, title1.c_str(), title2.c_str(), seqLength1, seqLength2, fpout);
                WriteAlignmentNeedle(title1.c_str(), title2.c_str(), alignXstr, alignYstr, alignRel, iAlign,nchar, fpout);
            } else if (outformat == 1){
                WriteAlignmentInFasta( title1, title2, alignXstr, gapopen, gapextend, iAlign, seqLength1, &alignFactor, fpout  );
                WriteAlignmentInFasta( title2, title1, alignYstr, gapopen, gapextend, iAlign, seqLength2, &alignFactor, fpout  );
            }
            if (fpTable != NULL) {
                WriteAlignmentTableLine(title1, title2, iAlign, seqLength1, seqLength2, alignFactor, fpTable);
            }
        }
        fclose(fpseq2);
    }
    return 0;
}/*}}}*/
int GenerateSequentialSelection(int &i, int &j, int sizeI, int sizeJ)/*{{{*/
{ /*generate i, j pair like a two layer for loop*/
    if (j < sizeJ -1  ){
        j ++;
    }else {
        i ++;
        j = i+1;
    }
    if (j < sizeJ && i < sizeI){
        return 0;
    }else {
        return -1;
    }
}/*}}}*/
int GenerateRandomSelection(int &i, int &j, int sizeI, int sizeJ, set<string>&selectedSet)/*{{{*/
{
    int max_iteration = 20;
    int cnt_iteration = 0;
    while(1){
        cnt_iteration ++;
        i = int (sizeI * double( rand())/( RAND_MAX+1.0));
        j = int (sizeJ * double( rand())/( RAND_MAX+1.0));
        if (i != j){
            string selectedpair = int2string (min(i,j)) + string("-") + int2string(max(i,j));
            if (selectedSet.find(selectedpair) == selectedSet.end()){ /*it is a new pair*/
                /*selectedSet.insert(selectedpair);*/ /*deleted 2011-10-31,
                                                        this may cause memory
                                                        overflow when the
                                                        comparison number is
                                                        big, e.g. 200000 x 200000*/
                return 0;
            }
        }
        if (cnt_iteration > max_iteration){
            break;
        }
    }
    return -1; /*after max_iteration, still no new random pairs are selected,
                 it means they are almost full, randomness will be bad*/
}/*}}}*/

int main(int argc, char *argv[])/*{{{*/
{	
    if(0){
        for(int i=0;i<numSeqIDTClass;i++){
            printf("%d: %g - %g\n",i, binSeqIDTClass[2*i], binSeqIDTClass[2*i+1]);
        }
    }

    if (argc < 2) {
        fprintf(stdout,"too few arguments\n");
        PrintHelp();
        return  1;
    }
    
    string programname = "my_needle";
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    string seqDatabase = "";/*indexed sequence database*/
    vector <string> positionalArgList; /*should be seqfile1 seqfile2*/

    string  outfile = "";
    string outTableFile = "";
    //GetDataDir(datadir);
    //sprintf(matrixFile, "%s/%s/%s", datadir, "matrices", "BLOSUM62");


    vector <Pair> pairList;
    pairList.clear();

    string pairListFile = "";
    string logFile = "";
    
    string title1 = "";
    string title2 = "";

    int i = 0 ;
    i = 1;
    bool isNonOptionArg = false;
    while(i < argc)/*{{{*/ {
        if (argv[i][0] == '-' && !isNonOptionArg) {
            if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                PrintHelp();
                return 0;
            } else if(string(argv[i])== "-m" || string(argv[i])== "--m" || string(argv[i]) == "-outformat" || string(argv[i]) == "--outformat") {
                outformat = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i], "-mode") == 0 || strcmp(argv[i], "--mode")== 0) {
                pairlistmode = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i], "-seqdb") == 0 || strcmp(argv[i], "--seqdb")== 0) {
                seqDatabase=argv[i+1];
                i += 2;
            } else if(strcmp(argv[i], "-G") == 0) {
                gapopen = atof(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i], "-E") == 0) {
                gapextend = atof(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i], "-endweight") == 0) {
                endweight = ajTrue;
                i += 1;
            } else if(strcmp(argv[i], "-eg") == 0) {
                endgapopen = atof(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i], "-ee") == 0) {
                endgapextend = atof(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i], "--nchar") == 0) {
                nchar = atoi(argv[i+1]);
                i += 2;
            } else if(string(argv[i])== "-type" || string(argv[i])== "--type") {
                alnType = atoi(argv[i+1]);
                if (alnType < 0 ||  alnType > 2) {
                    fprintf(stderr,"Wrong alignment type:'%d'! The alignment type should be 0 or 1 or 2\n", alnType);
                }
                i += 2;
            } else if(strcmp(argv[i], "-matrix") == 0 || strcmp(argv[i], "--matrix") == 0) {
                matrixFile=argv[i+1];
                i += 2;
            } else if(strcmp(argv[i], "-l" ) == 0 || strcmp(argv[i], "-list") == 0 || strcmp(argv[i], "--list") == 0) {
                pairListFile = argv[i+1];
                i += 2;
            } else if(strcmp(argv[i], "-title1") == 0 || strcmp(argv[i], "--title1") == 0) {
                title1 = argv[i+1];
                i += 2;
            } else if(strcmp(argv[i], "-title2") == 0 || strcmp(argv[i], "--title1") == 0) {
                title2 = argv[i+1];
                i += 2;
            } else if(strcmp(argv[i], "-debug") == 0 || strcmp(argv[i], "--debug") == 0) {
                logFile = argv[i+1];
                isPrintTraceMatrix = true;
                show = true;
                i += 2;
            } else if(string(argv[i])== "-threshold-seqdb-size" || 
                    string(argv[i]) == "--threshold-seqdb-size" )  {
                THRESHOLD_DB_FILESIZE = atoi(argv[i+1]);
                i += 2;
            } else if(string(argv[i])== "-max-align-size" || 
                    string(argv[i]) == "--max-align-size" )  {
                MAX_ALIGN_SIZE = atoi(argv[i+1]);
                i += 2;
            } else if(string(argv[i])== "-selpair" || 
                    string(argv[i]) == "--selpair" )  {
                method_select_pair = atoi(argv[i+1]);
                i += 2;
            } else if(string(argv[i])== "-wordsize" || 
                    string(argv[i]) == "--wordsize" )  {
                wordsize = atoi(argv[i+1]);
                i += 2;
            } else if(string(argv[i])== "-pkmer" || 
                    string(argv[i]) == "--pkmer" ||
                    string(argv[i]) == "-printkmerscore" ||
                    string(argv[i]) == "--printkmerscore" 
                    )  {
                if (toupper(argv[i+1][0]) == 'Y'){
                    isPrintKMerBitScore = true;
                }else{
                    isPrintKMerBitScore = false;
                }
                i += 2;
            } else if(string(argv[i])== "-o" || string(argv[i]) == "--o") {
                outfile = argv[i+1];
                i += 2;
            } else if(string(argv[i])== "-table" || string(argv[i]) == "--table") {
                outTableFile = argv[i+1];
                i += 2;
            } else if(strcmp(argv[i], "-nogoarray") == 0 || strcmp(argv[i], "--nogoarray") == 0) {
                noGapOpenArray = true;
                i += 1;
            } else if (strcmp(argv[i], "--") == 0)/*next item is non option argument*/ {
                isNonOptionArg = true;
                i ++;
                continue;
            } else {
                fprintf(stderr,"Error! Invalid argument '%s'\n", argv[i]);
                return -1;
            }
        } else {
            positionalArgList.push_back(argv[i]);
            if (positionalArgList.size() > 2){
                fprintf(stderr,"Error! argument %s is not valid. only two positional arguments can be supplied\n", argv[i]);
                return -1;
            }
            i ++;
        }
    }/*}}}*/

    /*Conver size from MB to Byte*/
    MAX_ALIGN_SIZE *= 1024*1024;  
    THRESHOLD_DB_FILESIZE *= 1024*1024;

    if (method_select_pair == 0){
        isRandSelection = false;
    }else{
        isRandSelection = true;
    }

    int numPair = 0;
    string multiseqFastaFile = "";
    if (pairlistmode == 0 || pairlistmode == 1 ){
        if (positionalArgList.size() == 2){
            Pair tmppair;
            tmppair.mem1 = positionalArgList[0];
            tmppair.mem2 = positionalArgList[1];
            pairList.push_back(tmppair);
        } else if (positionalArgList.size() != 0) {
            fprintf(stderr, "Exactly two positional arguments are accepted, but %d has been supplied, when mode is 1 or 2. Exit.\n", (int)positionalArgList.size());
            return 1;
        }
        if (pairListFile != ""){
            ReadInPairList(pairListFile.c_str(), pairList);
        }
        numPair = pairList.size();
        if (numPair <= 0){
            fprintf(stderr,"No pair set. Exit.\n");
            return 1;
        }
    }else if (pairlistmode == 2){
        if (positionalArgList.size () == 1){
            multiseqFastaFile = positionalArgList[0];
        }else{
            cerr << "Only one (but " << positionalArgList.size() 
                << " were supplied) positional argument should be set when mode is 2. Exit."
                << endl;
            return 1;
        }
    } else {
        cerr << "mode " << pairlistmode << " not implemented yet. Exit." << endl;
        return 1;
    }

    if (!endweight) {
        endgapopen = 0.0;
        endgapextend = 0.0;
    }

    FILE *fpout = stdout;
    if(outfile !=  "") {
         fpout = fopen(outfile.c_str(),"w");
         if (fpout == NULL){
             cerr << "Failed to write to file " << outfile 
                 << ". Reset output to STDOUT." << endl;
             fpout = stdout;
         }
    }

    FILE *fpTable = NULL;
    if (outTableFile != ""){
        fpTable = fopen(outTableFile.c_str(),"w");
        if (fpTable == NULL){
            cerr << "Failed to write to table file " << outTableFile << endl;
        }
    }

    FILE *fpLog = NULL;
    if(logFile !=  "" ) {
        fpLog = fopen(logFile.c_str(), "w");
        if (fpLog == NULL) {
            fprintf(stderr,"Failed to write to log file  %s.\n", logFile.c_str());
        }
    }

    int maxSizeSubMatrix = max(NUM_BLOSUM,max(NUM_NUC, NUM_SHAG));
    Array2D <float> SM_2darray(maxSizeSubMatrix+1,maxSizeSubMatrix+1);
    SM_2darray.Init(0);
    float **SM = SM_2darray.array2D;
    if (GetSubstitutionMatrix(matrixFile, SM, alphabet, alnType) == -1) {
        cerr << "Fatal! Get substitution matrix failed. Exit." << endl;
        exit(1);
    }

    if(outformat != 0 && outformat != 1) {
        fprintf(stderr,"Error! outformat = %d not implemented. Reset to 0.", outformat);
        outformat = 0;
    }  
    if (outformat==0){
        /*write the overal header for needle output*/
        PrintNeedleAlignmentHeader(programname, timeinfo,outfile, argc, argv,fpout); 
    } 
    if (fpTable != NULL){
        fprintf(fpTable,"%s", explanation_table_format);
        fprintf(fpTable,"#%-15s %-15s %6s %6s %9s %6s %6s %9s %6s %6s %6s %6s %6s\n", "Seq1","Seq2", "IDT0", "SIM0", "AlnLength", "Len1","Len2", "Score","N_IDT", "N_SIM", "N_GAP", "IDT1", "IDT2");
    } 

    if (pairlistmode == 0){   /*pair of files */
        for (i=0 ; i < numPair ; i++){ 
            AlignTwoSeqFile(pairList[i].mem1, pairList[i].mem2, title1, title2, SM, fpout, fpLog, fpTable);
        }
    }else if (pairlistmode == 1){
        if (seqDatabase == ""){
            cerr << "seqDatabase not set. It's required when pairlistmode = 1. Exit." << endl;
            return -1;
        }
        map <string, dbindex> dbindexmap;
        int maxDBIndexNumber = 0;
        /*=== Begin Read DB index */
#ifdef DEBUG_TIME
        clock_t start,finish;
        double duration;
        start=clock();
#endif
        if (ReadDatabaseIndex(seqDatabase, dbindexmap, maxDBIndexNumber) == -1 ){
            cerr << "Read Database Index failed for "<<seqDatabase << ". Exit." << endl;
            return -1;
        }
#ifdef DEBUG_TIME
        finish=clock();
        duration = double(finish-start)  /double(CLOCKS_PER_SEC);
        cout << "Read Database Index with " << dbindexmap.size()<< " records costs " << duration << " seconds." << endl;
#endif
        vector <FILE*> fpList;
        GetDBFPList( fpList, seqDatabase, maxDBIndexNumber);
        /*=== End Read DB index */
#ifdef DEBUG
        cout << "Read index finished" << endl;
#endif
        int dbfilesize = GetFileSize(seqDatabase.c_str());
        map <string, string> idseqmap;
        if (dbfilesize < THRESHOLD_DB_FILESIZE){
#ifdef DEBUG_TIME
            start=clock();
#endif
            /*Read the whole seqDatabase into memory*/
            GetIDSeqMap(dbindexmap,fpList, idseqmap);
#ifdef DEBUG_TIME
            finish=clock();
            duration = double(finish-start)  /double(CLOCKS_PER_SEC);
            cout << "GetIDSeqMap of "<<idseqmap.size()<<" sequences costs " <<
                duration << " seconds." << endl;
#endif
        }
#ifdef DEBUG_TIME
        start=clock();
#endif
        AlignAllPair(pairList, dbindexmap, fpList, idseqmap, dbfilesize, SM, fpout, fpLog, fpTable);
        for(int i=0;i <= maxDBIndexNumber; i++){
            fclose(fpList[i]);
        }
#ifdef DEBUG_TIME
        finish=clock();
        duration = double(finish-start)  /double(CLOCKS_PER_SEC);
        cout << "Align "<<pairList.size()<<" pairs of sequences costs " <<
            duration << " seconds." << endl;
#endif
    } else if (pairlistmode == 2) {
        if (multiseqFastaFile == ""){
            cerr << "Fasta file not set for mode " << pairlistmode << ". Exit" << endl;
            return 1;
        }
        int numSeq = 0;
        vector <string> idList;
        vector <string> seqList;
        if (isShowProgress){
            cout << "Reading fasta file " << multiseqFastaFile  << endl ;
        }
#ifdef DEBUG_TIME
        clock_t start, finish;
        double duration;
        start = clock();
#endif
        numSeq = ReadFasta(multiseqFastaFile, idList, seqList, seq_type) ;
#ifdef DEBUG_TIME
        finish = clock();
        duration = double(finish-start)  /double(CLOCKS_PER_SEC);
        cout << "Read "<< numSeq <<" sequences costs " <<
            duration << " seconds." << endl;
#endif
        if (isShowProgress){
            cout << numSeq << " sequences have been read in to memory." << endl ;
        }
        if (numSeq < 2){
            cerr << "Too few sequences (" << numSeq << ") read from file " << multiseqFastaFile
                << ". Exit." << endl;
            return 1;
        }
        long long cntAllOutputAln = 0;
        long long cntAligned = 0;
        long long numTotalPair =  ((long long)numSeq*(numSeq-1)/2) ;
        vector <int> cntSeqIDTClass;
        for ( int kk = 0; kk < numSeqIDTClass; kk ++ ){
            cntSeqIDTClass.push_back(0);
        }
        if (isShowProgress){
            cout << "Start all-to-all pairwise alignment (" <<  ((long long)numSeq*(numSeq-1)/2) << " alignments)..." << endl ;
        }
        int i = 0;
        int j = 0;
        set <string> selectedSet;
        if(isRandSelection){
            srand(time(NULL));
        }

        long long cntPairRejectedByKMerCompare=0;
        int sizeTable = int(pow(double(23), double(wordsize)));
        int *freqTable1 = new int [sizeTable+1];
        int *freqTable2 = new int [sizeTable+1];
#ifdef DEBUG_TIME
        start = clock();
#endif
        float minimalSeqIdentity = 0.0;
        float threshold_kmerScore = -99999;
        bool isRunKMerComparison = false;
        long long cntPair = 0;
        while (1) {
            if (!isRandSelection){
                if (GenerateSequentialSelection(i,j, numSeq, numSeq) != 0){
                    break;
                }
            }else{
                if (GenerateRandomSelection(i,j, numSeq, numSeq, selectedSet) != 0){
                    break;
                }
            }
            cntPair ++;
            if (cntPair > numTotalPair){
                fprintf(stdout,"cntPair = %lld. Break the loop.\n",cntPair);
                break;
            }
#ifdef DEBUG_SELECTION
            fprintf(stdout, "pair %lld : %d - %d\n", cntPair, i,j);
            continue;
#endif

#ifdef DEBUG_KMER_BITSCORE
            minimalSeqIdentity = 100.0; /*DEBUG*/
#endif
            if (minimalSeqIdentity >= 10.0) {
                isRunKMerComparison = true;
            }
            
            float kmerBitScore = 0.0;
            if (isCheckSeqIDTClass){
                minimalSeqIdentity = GetMinimalSeqIdentity(cntSeqIDTClass,
                        maxSeqIDTClass, binSeqIDTClass, numSeqIDTClass);
                threshold_kmerScore =
                    GetThresholdKMerScore(minimalSeqIdentity);
            }
            //printf("minimalSeqIdentity = %g threshold_kmerScore = %g\n", minimalSeqIdentity, threshold_kmerScore);
            if (isRunKMerComparison){
                /*do faster comparison*/
                //kmerBitScore = KMerVectorPairwiseComparison_c(seqList[i], seqList[j], idList[i], idList[j], wordsize, fpLog);
                //kmerBitScore = KMerVectorPairwiseComparison(seqList[i], seqList[j], idList[i], idList[j], wordsize, fpLog);
                kmerBitScore = KMerVectorPairwiseComparison_table(seqList[i], seqList[j], idList[i], idList[j], freqTable1, freqTable2, sizeTable, wordsize, fpLog);
                if (kmerBitScore < threshold_kmerScore){
                    cntPairRejectedByKMerCompare ++;
                    continue;
                }
            }
            //continue;
            int iAlign = 0;/*alignment length*/
            int seqLength1 = seqList[i].size();
            int seqLength2 = seqList[j].size();
            Array1D <char> alignXstr_1darray(seqLength1+seqLength2+1);
            Array1D <char> alignYstr_1darray(seqLength1+seqLength2+1);
            Array1D <int>  alignRel_1darray(seqLength1+seqLength2+1);
            char *alignXstr = alignXstr_1darray.array1D;
            char *alignYstr = alignYstr_1darray.array1D;
            int  *alignRel  = alignRel_1darray.array1D;
            AlignFactor alignFactor;
            InitAlignFactor(&alignFactor);
            iAlign = Alignment_MyNeedle(seqList[i].c_str(), seqList[j].c_str(), alphabet, seqList[i].size(), seqList[j].size(), idList[i].c_str(), idList[j].c_str(), alignXstr,alignYstr,alignRel, alignFactor, gapopen, gapextend, endgapopen, endgapextend, SM, seq_type, fpLog);
            cntAligned += 1;

            if(isRunKMerComparison && isPrintKMerBitScore) {
            fprintf(stdout,"%-12s - %12s vector bitScore = %10.3g seqidt = %10.3g len1 = %d len2 = %d iAlign = %d\n", idList[i].c_str(),idList[j].c_str(),kmerBitScore, alignFactor.identity*100, seqLength1, seqLength2, iAlign);
            }

            if (isShowProgress) {
                if (cntAligned%1000 == 1){ /*write a newline will flush the buffer*/
                    fprintf(stdout,"Processed alignments %lld, numpair rejected by KMer compare %lld.\n", cntAligned, cntPairRejectedByKMerCompare);
                    if (isCheckSeqIDTClass) {
                        for ( int i = 0 ; i < numSeqIDTClass; i ++){
                            fprintf(stdout,"\tcntSeqIDTClass[%d] = %d\n",i, cntSeqIDTClass[i]);
                        }
                        fprintf(stdout,"\tTotal output alignments = %lld\n",cntAllOutputAln);
                    }
#ifdef DEBUG_TIME
                    finish = clock();
                    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
                    cout << "Time used for comparison: " <<
                        duration << " seconds." << endl;
#endif
                    fflush(stdout);
                }
            }
            int seqidtclass = 0;
            if (isCheckSeqIDTClass){
                seqidtclass = GetSeqIDTClass(alignFactor.identity*100, binSeqIDTClass, numSeqIDTClass);
#ifdef DEBUG_SEQIDTCLASS
                fprintf(stdout,"seqidt %.3g seqidtclass = %d\n", alignFactor.identity*100, seqidtclass);
#endif
            }
            if ( !isCheckSeqIDTClass || (seqidtclass != numSeqIDTClass && cntSeqIDTClass[seqidtclass] < maxSeqIDTClass[seqidtclass]) ) {
                if (outformat == 0){
                    WriteAlignmentHeaderNeedle(gapopen, gapextend, iAlign, &alignFactor, idList[i].c_str(), idList[j].c_str(), seqLength1, seqLength2, fpout);
                    WriteAlignmentNeedle(idList[i].c_str(), idList[j].c_str(), alignXstr, alignYstr, alignRel, iAlign,nchar, fpout);
                }else if (outformat==1){
                    WriteAlignmentInFasta( idList[i], idList[j], alignXstr, gapopen, gapextend, iAlign, seqLength1, &alignFactor, fpout  );
                    WriteAlignmentInFasta( idList[j], idList[i], alignYstr, gapopen, gapextend, iAlign, seqLength2, &alignFactor, fpout  );
                }
                if (fpTable != NULL){
                    WriteAlignmentTableLine(idList[i], idList[j], iAlign, seqLength1, seqLength2, alignFactor, fpTable);
                }
                cntSeqIDTClass[seqidtclass] +=1;
                cntAllOutputAln += 1;
                string selectedpair = int2string (min(i,j)) + string("-") + int2string(max(i,j));
                selectedSet.insert(selectedpair);/*added selected pair when only it actually output*/
            }
        }
        delete [] freqTable1;
        delete [] freqTable2;
#ifdef DEBUG_TIME
        finish = clock();
        duration = double(finish-start)  /double(CLOCKS_PER_SEC);
        cout << "Compare "<< cntPair <<" pairs costs " <<
            duration << " seconds." << endl;
#endif
        if (isShowProgress) {
            fprintf(stdout,"Finished.\n");
        }
        fprintf(stdout,"Total processed alignments %lld, numpair rejected by KMer compare %lld.\n", cntAligned, cntPairRejectedByKMerCompare);
        if (isCheckSeqIDTClass) {
            for ( int i = 0 ; i < numSeqIDTClass; i ++){
                fprintf(stdout,"\tcntSeqIDTClass[%d] = %d\n",i, cntSeqIDTClass[i]);
            }
            fprintf(stdout,"\tTotal output alignments = %lld\n",cntAllOutputAln);
        }
    } else{
        fprintf(stderr,"pairlistmode %d not implemented yet. Exit", pairlistmode);
        return 1;
    }
    
    /*write the end of the file*/
    if (outformat == 0){
        fprintf(fpout,"#---------------------------------------\n");
        fprintf(fpout,"#---------------------------------------\n");
    }

    if(fpout != stdout && fpout != NULL) {
        fclose(fpout);
    }
    if (fpLog != NULL && fpLog != stdout) {
        fclose(fpLog);
    }
    if (fpTable != NULL && fpTable != stdout) {
        fclose(fpTable);
    }
	return 0;
}
/*}}}*/
