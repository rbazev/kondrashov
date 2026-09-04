# Kondrashov

Attempt to repeat the following study:

[A. S. Kondrashov, S. Sunyaev, and F. A. Kondrashov. Dobzhansky–Muller
incompatibilities in protein evolution. *Proc. Natl. Acad. Sci. U. S. A.*
99(23): 14878–83, 2002.](https://www.pnas.org/doi/10.1073/pnas.232565499)







The analysis has the following components:

1. [orthologs.ipynb](orthologs.ipynb) collect primate orthologs of each human gene from NCBI, remove redundant sequences, save them to a `*.fasta` file.  

2. [Clinvar_pathogen_Mary.ipynb](Clinvar_pathogen_Mary.ipynb) extract pathogenic variants from ClinVar data, save as .csv file in 'pathogenic folder'

3. [find_closest.ipynb](find_closest_Mary.ipynb) Identify variants that are most similar to reference human protein within the various species

4. [align.ipynb](align.ipynb) align orthologous sequences in each `*.fasta` file and save it to a `*.aln` file.

5. [find_CPDs.ipynb](find_CPDs.ipynb) find putative CPDs.

6. [clean_up.ipynb](clean_up.ipynb) establish connections between transcripts in ClinVar data and proteins in ortholog collections.