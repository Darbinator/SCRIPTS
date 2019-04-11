import scipy.stats as stats
import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import re

parser = argparse.ArgumentParser(description='test of overrepresentation of chromatin state in a list of specific genes\n\nExemple: python3 fisher_chrState.py 3 DEG/Col_ColHS/id.txt', formatter_class=RawTextHelpFormatter)
parser.add_argument("state",action='store', type=int, help="integer for a chromatin state (1-36) \n\n\
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
parser.add_argument("--pvalue",action='store', help="p-value associate to the exact fisher test for overepresentation of chromatin state (default=0.01) ", default=0.01)

args = parser.parse_args()


len_genes_background = 27416


if os.path.isfile(args.list):
	with open(args.list, 'r') as fp:
	            args.list_read = fp.read()

	            regex_TAIRID = re.compile("AT[A-Z0-9]{1}G[0-9]{5}")

	            genes_input = regex_TAIRID.findall(args.list_read)


else:
	print("Please give a correct list of genes file")
	exit(0)
	
if args.state not in range(1,37):
	print("Please give a correct chromatin state number (between 1 and 36)")



else:
	if not os.path.exists("At_genes_S"+str(args.state)):
		os.system("wget http://systemsbiology.cau.edu.cn/chromstates/download/At_genes_S"+str(args.state))

	with open("At_genes_S"+str(args.state),"r") as file:

		target_background = []
		for lignes in file.readlines():
			target_background.append(lignes.split()[1].rstrip())


target_input = []
for code in target_background:
	if code in genes_input:
		target_input.append(code)




oddsratio, pvalue = stats.fisher_exact([[len_genes_background, len(target_background)], [len(genes_input), len(target_input)]])

if pvalue>args.pvalue:
	print(args.state,pvalue,"not sign")

elif oddsratio < 1:
	print(args.state,pvalue,"under")

else:
	print(args.state,pvalue,"over")

#os.path.dirname(args.list)






