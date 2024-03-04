# korthoCPP
korthoCPP calculates pairwise kmer jaccard distance between two peptide fasta files

## Download and compile

```
git clone https://github.com/kullrich/korthoCPP.git
cd korthoCPP && make
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
korthoCPP -q Pseudomonas_aeruginosa_PAO1_107.faa -t Pseudomonas_fluorescens_SBW25_116.faa -p 2
head output.txt
```
```
