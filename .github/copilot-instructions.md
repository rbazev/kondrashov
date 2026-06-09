# Copilot Instructions for Kondrashov Project

## Project Overview
This project replicates the 2002 Kondrashov et al. study analyzing protein evolution for Dobzhansky–Muller incompatibilities across 32 disease-associated human genes. It integrates primate sequence data from NCBI, functional variant information from ClinVar, and sequence alignment analysis.

## Architecture & Data Flow

### Five-Stage Analysis Pipeline
The project uses **Jupyter notebooks** for interactive analysis, each stage outputs data for the next:

1. **orthologs.ipynb** → Retrieves primate orthologs from NCBI; outputs to `fasta/{gene}.fasta`
2. **align.ipynb** → Aligns sequences using `muscle` command; outputs to `aln/{gene}.aln`
3. **pathogen.ipynb** → Extracts pathogenic variants from ClinVar `clinvar/` folder
4. **clean_up.ipynb** → Maps ClinVar transcripts to ortholog protein sequences
5. **find_CPDs.ipynb** → Identifies putative Dobzhansky–Muller incompatibilities (CPDs)

### Key Modules
- **kondrashov.py** (~914 lines): Core utilities split into 5 sections (see comments):
  - **RefSeq Transcripts**: `get_transcript_ids()`, `get_protein_ids()`, sequence retrieval
  - **Sequence Processing**: CDS extraction, mutation analysis, codon translation
  - **ClinVar Variants**: Variant fetching, parsing phenotypic effects
  - **Alignment Analysis**: `get_alleles()`, allele frequency checks, species comparison
  - **Mutation Classification**: Synonym/nonsynonym/nonsense categorization

### External Dependencies
- **NCBI Entrez API** (requires credentials):
  - Email/API key configured in `kondrashov.py` (lines 2-3)
  - Used by all transcript/protein/variant lookup functions
- **BioPython**: Sequence I/O (SeqIO, AlignIO), translation, parsing
- **muscle**: External MSA tool required by align.ipynb (must be installed)
- **pandas, xmltodict**: Data handling and XML parsing

## Code Patterns & Conventions

### Naming & Structure
- **Gene list**: 32 loci hardcoded as `loci` tuple in kondrashov.py (ABCD1, ALPL, AR, ...)
- **Amino acid mapping**: `one_letter` dict (Ala→A); `nucl` tuple defines DNA bases
- **Protein change format**: ClinVar entries parsed as "p.XyzNNMut" (3-letter codes)
- **Transcript IDs**: NM_* (nucleotide) vs NP_* (protein); functions filter by type

### File I/O Patterns
- FASTA/alignment files named by gene: `fasta/{gene}.fasta`, `aln/{gene}.aln`
- ClinVar data in `clinvar/` and `pathogenic/` subdirectories (checked by pathogen.ipynb)
- All file paths are **relative to project root** when run in notebooks

### API Usage Caution
- `get_transcript_ids()`, `get_protein_ids()`: Filter results (return lists, not raw search)
- `get_all_changes()`: Reads local CSV/TXT files, not live API—depends on pathogen.ipynb preprocessing
- Error handling: Most functions catch exceptions silently (e.g., variant parsing)

## Common Development Tasks

### Adding Analysis for New Gene
1. Add gene to `loci` tuple in kondrashov.py
2. Rerun orthologs.ipynb (fetches NCBI data → fasta/)
3. Rerun align.ipynb (generates → aln/)
4. Rest of pipeline auto-includes new gene

### Extending Mutation Analysis
- Edit `get_codon_mutations()` or `get_mutation_counts()` in kondrashov.py
- These are called from notebooks only; test with small gene subset first
- Remember: nonsense mutations (→ stop codon) tracked separately from missense

### Debugging Variant Parsing
- `get_change_from_variant()` assumes format "p.XyzNNMut"—fails silently on malformed entries
- `get_variant_details()` parses ClinVar XML structure; different record types handled separately
- Check actual ClinVar data format if functions return incomplete results

## Notebook Execution Notes
- Notebooks import `kondrashov` module with `run kondrashov` magic command
- Execution order matters: later notebooks depend on earlier outputs (fasta, aln files)
- orthologs.ipynb has pre-computed outputs (some cells show stdout without re-running)
- NCBI queries in orthologs.ipynb can be slow; already-fetched data cached locally

## Key Functions for Common Tasks

| Function | Purpose | Input | Output |
|----------|---------|-------|--------|
| `get_alleles()` | Find all species alleles at a site | gene, protein ID, amino acid site | bool (variant exists), site, human AA, counts |
| `get_all_changes()` | Parse all variants for a transcript | gene, transcript ID, pathogenic_only flag | lists: sites, wilds, muts, parsing errors |
| `get_mutation_counts()` | Classify all single mutations | DNA sequence | (synonymous, nonsynonymous, nonsense) counts |
| `global_compare_to_human()` | Identify differences between aligned sequences | gene, human ID, species ID | dict of differences by position |
