# Kondrashov

Attempt to repeat the following study:

Original study
A. S. Kondrashov, S. Sunyaev, and F. A. Kondrashov. Dobzhansky–Muller incompatibilities in protein evolution. *Proc. Natl. Acad. Sci. U. S. A.* 99(23): 14878–83, 2002.

The analysis has the following components:

1. [orthologs.ipynb](orthologs.ipynb) collect primate orthologs of each human gene from NCBI, remove redundant sequences, save them to a `*.fasta` file.  

1. [align.ipynb](align.ipynb) align orthologous sequences in each `*.fasta` file and save it to a `*.aln` file.

1. [pathogen.ipynb](pathogen.ipynb) extract pathogenic variants from ClinVar data. 

1. [clean_up.ipynb](clean_up.ipynb) establish connections between transcripts in ClinVar data and proteins in ortholog collections.

1. [find_CPDs.ipynb](find_CPDs.ipynb) find putative CPDs.