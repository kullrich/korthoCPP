# korthoCPP
korthoCPP calculates pairwise kmer jaccard distance between two peptide fasta files

## Download and compile

```
git clone https://github.com/kullrich/korthoCPP.git
cd korthoCPP && make clean && make
```

## Usage

```
Usage: korthoCPP [options] -q <query.fasta> -t <target.fasta>
Options:
  -q FILE    query peptide fasta file
  -t FILE    target peptide fasta file
  -o FILE    output file to write jaccard distances
  -k INT     kmer length [default: 6]
  -m DOUBLE  min jaccard distance to report pair [default: 0.01]
  -s DOUBLE  sparse threshold to switch search strategy [defualt: 0.1]
  -n INT     number of kmers to check for sparse [default: 20]
  -p INT     number of threads [default: 1]
  -d         debug
```

## Compare two peptide fasta files and calcualte jaccard distances

```
korthoCPP -q query.fasta -t target.fasta
```

### use multiple threads

```
korthoCPP -q query.fasta -t target.fasta -p 2
```

### use different min jaccard

```
korthoCPP -q query.fasta -t target.fasta -m 0.02
```

## Example

```
wget https://www.pseudomonas.com/downloads/pseudomonas/pgd_r_22_1/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107.faa.gz
wget https://www.pseudomonas.com/downloads/pseudomonas/pgd_r_22_1/Pseudomonas_fluorescens_SBW25_116/Pseudomonas_fluorescens_SBW25_116.faa.gz
gunzip Pseudomonas_aeruginosa_PAO1_107.faa.gz
gunzip Pseudomonas_fluorescens_SBW25_116.faa.gz
./korthoCPP -q Pseudomonas_aeruginosa_PAO1_107.faa -t Pseudomonas_fluorescens_SBW25_116.faa -p 4 -d
```

```
Time taken: kmer extraction 7059 milliseconds
number of sequences Q: 5586
number of sequences T: 5861
Time taken: QkmerMap creation 7091 milliseconds
number of kmers QkmerMap: 1612859
Time taken: TkmerMap creation 8469 milliseconds
number of kmers TkmerMap: 1734068
number of sparse candidate pairs: 345
sparse value: 1.05377e-05 < sparse threshold: 0.1 >>> search strategy one vs one
Time taken: check sparse threshold 1 milliseconds
Time taken: candidate pairs creation 1423 milliseconds
number of candidate pairs: 409274
Time taken: distance calculation 44714 milliseconds
```

```
head output.txt
###
qname	tname	jaccard	mash	ani	sumdist
PA0034	PFLU1157	0.0303797	0.471793	0.528207	0.941032
PA0032	PFLU0028	0.087344	0.304749	0.695251	0.839344
PA0029	PFLU0025	0.057554	0.369641	0.630359	0.891156
PA0030	PFLU0026	0.0470383	0.401602	0.598398	0.91015
PA0031	PFLU0027	0.215854	0.172576	0.827424	0.644935
PA0035	PFLU0035	0.272947	0.141111	0.858889	0.571157
PA0036	PFLU0036	0.473297	0.0737314	0.926269	0.3575
PA0037	PFLU0037	0.110476	0.269099	0.730901	0.801029
PA0038	PFLU0039	0.147826	0.226074	0.773926	0.742424
```
