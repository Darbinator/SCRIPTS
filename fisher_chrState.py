import scipy.stats as stats
import sys


import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='test of overrepresentation of chromatin state in a list of specific genes', formatter_class=RawTextHelpFormatter)
parser.add_argument("state",action='store', help="integer for a chromatin state (1-36) \n\n\
list of chromatin states possible:\n\
state	Preferential Epigenetics Marks		Preferential location\n\
1			H3.3						3'UTR\n\
2			H3.3 histone acetylation,H3K4me2,H2A.Z		CDS,3'UTR\n\
3			H3K4me1,H3.3,H3.1				CDS\n\
4			H3K4me1,H3.3					CDS,intron\n\
5			H3K4me1,H3K36me3,H3.3,H3.1			CDS\n\
6			H3K4me1,H3K36me3				intron\n\
7			H3K4me1,H3K36me3,H3K4me2			CDS,intron\n\
8			H3K4me1,H3K4me2,H2A.Z				CDS\n\
9			H3K4me1						intron\n\
10		H2A.Z						CDS,intron\n\
11		H3K27me3,H2A.Z,H3K4me2				CDS\n\
12		H3K27me3,H2A.Z					Promoter,CDS,intron,intergenic\n\
13		H3K27me3,H2A.Z					Promoter,intergenic\n\
14		H3K27me3					Promoter,intergenic\n\
15		H3K27me3,accessible DNA				Promoter,intergenic\n\
16		accessible DNA					Promoter,intergenic\n\
17		accessible DNA					Promoter\n\
18		accessible DNA					Promoter\n\
19		accessible DNA					Promoter,snRNA\n\
20		accessible DNA 					Promoter\n\
21		accessible DNA					Promoter\n\
22		histone acetylation,H3K4me2			coding gene,miRNA,snoRNA\n\
23		accessible DNA,H3K36ac,H3K56ac,H4K16ac,H3K4me3	Promoter,5'UTR\n\
24		accessible DNA,histone acetylation,H3K4me3	5'UTR,snoRNA\n\
25		histone acetylation,H3K4me3,H3K4me2,H2A.Z	intron\n\
26		histone acetylation,H3K4me3,H3K4me2,H2A.Z	CDS\n\
27		H3K4me2,histone acetylation,H3K4me3,H2A.Z	CDS\n\
28		H3K4me3,H3K4me2,H2A.Z				intron\n\
29		weak signal					intergenic\n\
30		rare signal					intergenic\n\
31		DNA methylation,H3K9me2,H3K27me3		intergenic,miRNA\n\
32		DNA methylation,H3K9me2				intergenic,TE\n\
33		H3K9me2,DNA methylation				TE\n\
34		H3K9me2,DNA methylation,H3K27me1		TE,miRNA\n\
35		H3K9me2,DNA methylation,H2A.X			intergenic,pericentromere\n\
36		CENH3,H3K9me2,DNA methylation,accessible DNA	rRNA,tRNA,centromere\n\n")
parser.add_argument("list", action='store', help="list of genes to analyse overrepresentation of epigenetic mark")

args = parser.parse_args()


with open(args.list, 'r') as fp:
            args.list = fp.read().strip(' ,\n')

print(args.list)


oddsratio, pvalue = stats.fisher_exact([[28775, 442], [835, 35]])

