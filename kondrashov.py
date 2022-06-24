from datetime import date
from io import StringIO


from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline


# How to get a key is described here: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
Entrez.email = 'razevedo@uh.edu'
Entrez.api_key = '41b789f82f616365ce6551f9105906c1cb08'


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xmltodict as x2d


########################
### GLOBAL VARIABLES ###
########################


# nucleotides in DNA
nucl = ('a', 'c', 'g', 't')


# 3-letter to 1-letter amino acid codes
one_letter={'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*'}


# genes in Kondrashov study
loci = ['ABCD1', 'ALPL', 'AR', 'ATP7B', 'BTK', 'CASR', 'CBS', 'CFTR', 'CYBB',
        'F7', 'F8', 'F9', 'G6PD', 'GALT', 'GBA', 'GJB1', 'HBB', 'HPRT1',
        'IL2RG', 'KCNH2', 'KCNQ1', 'L1CAM', 'LDLR', 'MPZ', 'MYH7', 'TYR',
        'PAH', 'PMM2', 'RHO', 'TP53', 'TTR', 'VWF']


##########################
### RefSeq transcripts ###
##########################


def get_transcript_ids(gene):
    '''
    Get RefSeq transcript sequence IDs from gene name.

    Collect all curated transcript variants (i.e., only NM_* accessions).

    Parameters
    ----------
    gene : str
        Gene name.

    Return
    ------
    list
        NCBI record IDs.
    '''
    search_term = "Homo sapiens[Organism] AND {0}[Gene Name] AND biomol_mrna[PROP] AND refseq[filter]".format(gene)
    handle = Entrez.esearch(db="nucleotide", retmax=100, term=search_term, idtype="acc")
    search = Entrez.read(handle)
    search_ids = []
    for i in search['IdList']:
        if i[:2] == 'NM':
            search_ids.append(i)
    handle.close()
    return search_ids


def get_protein_ids(gene):
    '''
    Get RefSeq transcript sequence IDs from gene name.

    Collect all curated transcript variants (i.e., only NM_* accessions).

    Parameters
    ----------
    gene : str
        Gene name.

    Return
    ------
    list
        NCBI record IDs.
    '''
    search_term = "Homo sapiens[Organism] AND {0}[Gene Name] AND refseq[filter]".format(gene)
    handle = Entrez.esearch(db="protein", retmax=100, term=search_term, idtype="acc")
    search = Entrez.read(handle)
    search_ids = []
    for i in search['IdList']:
        search_ids.append(i)
    handle.close()
    return search_ids


def get_transcript(transcript):
    '''
    Get transcript record from ID.

    Parameters
    ----------
    transcript : str
        NCBI record ID.

    Return
    ------
    str
        Transcript record.
    '''
    handle = Entrez.efetch(db="nucleotide", id=transcript, rettype="gb", retmode="xml")
    record = Entrez.read(handle)
    handle.close()
    return record[0]


def get_protein(protein):
    '''
    Get protein record from ID.

    Parameters
    ----------
    protein : str
        NCBI record ID.

    Return
    ------
    str
        Protein record.
    '''
    handle = Entrez.efetch(db="protein", id=protein, rettype="gb", retmode="xml")
    record = Entrez.read(handle)
    handle.close()
    return record[0]


def get_transcript_from_protein(protein):
    """
    Get transcript ID from protein record.

    Parameters
    ----------
    protein : str
        NCBI record ID.

    Returns
    -------
    str
        NCBI record ID.
    """
    rec = get_protein(protein)
    for feat in rec['GBSeq_feature-table']:
        featkey = feat['GBFeature_key']
        if featkey=='CDS':
            quals = feat['GBFeature_quals']
            for qual in quals:
                qualkey = qual['GBQualifier_name']
                if qualkey=='coded_by':
                    cds = qual['GBQualifier_value']
                    return cds.split(':')[0]


def get_transcript_from_variant(variant):
    """
    Get transcript ID from ClinVar record.

    Parameters
    ----------
    variant : str
        Entry in 'Name' column for a gene obtained from ClinVar.

    Returns
    -------
    str
        NCBI record ID.
    """
    transcriptgene = variant.split(':')[0]
    transcript = transcriptgene.split('(')[0]
    return transcript


def get_change_from_variant(variant):
    """
    Interpret protein change from ClinVar record.

    Parameters
    ----------
    variant : str
        Entry in 'Name' column for a gene obtained from ClinVar.

    Returns
    -------
    tuple
        site : int
            Amino acid site (starting from 1).
        wild : str
            Wild type amino acid.
        mut : str
            Mutant amino acid.
    """
    change = transcriptgene = variant.split(' (p.')[1][:-1]
    wild = change[:3]
    mut = change[-3:]
    site = int(change[3:-3])
    return site, one_letter[wild], one_letter[mut]


def get_all_changes(gene, transcript, pathogenic_only):
    '''
    Retrieve all changes for a given gene linked to a particular transcript.

    Parameters
    ----------
    gene : str
        Gene name.
    transcript : str
        Transcript NCBI ID.
    pathogenic_only : bool
        Whether to consider only pathogenic variants.

    Returns
    -------
    tuple
        sites : list
            Mutant sites.
        wilds : list
            Wild type amino acids.
        muts : list
            Mutant amino acids.
        errs : list
            Names of variants where no change was detected.
    '''
    sites = []
    wilds = []
    muts = []
    errs = []
    if pathogenic_only:
        vardata = pd.read_csv('pathogenic/' + gene + '_clinvar.csv')
    else:
        vardata = pd.read_table('clinvar/' + gene + '_clinvar.txt.txt', sep='\t')
    n = len(vardata)
    for i in range(n):
        name = vardata['Name'][i]
        tmp_transcript = get_transcript_from_variant(name)
        if tmp_transcript==transcript:
            try:
                site, wild, mut = get_change_from_variant(name)
                sites.append(site)
                wilds.append(wild)
                muts.append(mut)
            except:
                accession = vardata['Accession'][i]
                errs.append((name, accession))
    return sites, wilds, muts, errs


def get_transcripts_from_variants(gene, pathogenic_only):
    """
    Get transcript ID from all ClinVar records for a gene.

    Parameters
    ----------
    gene : str
        Gene name.
    pathogenic_only : bool
        Whether to consider only pathogenic variants.[clean_up.ipynb](clean_up.ipynb)

    Returns
    -------
    list
        NCBI record IDs.
    """
    if pathogenic_only:
        vardata = pd.read_csv('pathogenic/' + gene + '_clinvar.csv')
    else:
        vardata = pd.read_table('clinvar/' + gene + '_clinvar.txt.txt', sep='\t')
    vartranscripts = []
    names = vardata['Name'].tolist()
    for name in names:
        transcript = get_transcript_from_variant(name)
        vartranscripts.append(transcript)
    return list(set(vartranscripts))


def get_aln_positions(gene, protein):
    '''
    Get list of sites in alignment corresponding to sites in human sequence.

    The site j in position i in the output means that site i in the human protein is in position j in the alignment.

    Parameters
    ----------
    gene : str
        Gene name.
    protein : str
        Protein NCBI ID.

    Returns
    -------
    list
        Sites (starting at 0).
    '''
    for record in SeqIO.parse("fasta/{0}.fasta".format(gene), "fasta"):
        recordid = record.description.split(' ')[0]
        if recordid == protein:
            fasta = record
            break
    align = AlignIO.read('aln/{0}.aln'.format(gene), 'fasta')
    for record in align:
        if record.id == protein:
            aln = record
            break
    i = 0
    j = 0
    L = len(fasta)
    alni = []
    while i < L:
        if aln.seq[j] != '-':
            alni.append(j)
            i += 1
            j += 1
        else:
            j += 1
    return alni


def remove_element(x, element):
    '''
    Remove element from list.  If element is not in list do nothing.

    Parameters
    ----------
    x : list
        List.
    element : any
        List element to remove.
    '''
    try:
        i = x.index(element)
        del x[i]
    except ValueError:
        pass


def get_variable_sites(gene, protein, verbose):
    '''
    Identify sites that are variable among the sequences in the alignment.

    Exclude gaps.

    Parameters
    ----------
    gene : str
        Gene name.
    protein : str
        Protein NCBI ID.

    Returns
    -------
    dict
        Other alleles present at each site.
    '''
    variable = {}
    if verbose:
        print('Orig\tAln\tHuman\tOthers')
    for record in SeqIO.parse("fasta/{0}.fasta".format(gene), "fasta"):
        recordid = record.description.split(' ')[0]
        if recordid == protein:
            fasta = record
            break
    align = AlignIO.read('aln/{0}.aln'.format(gene), 'fasta')
    for record in align:
        if record.id == protein:
            aln = record
            break
    alni = get_aln_positions(gene, protein)
    L = len(fasta)
    for i in range(L):
        j = alni[i]
        counts = pd.Series(list(align[:, j])).value_counts()
        alleles = counts.index.tolist()
        remove_element(alleles, '-')
        n = len(alleles)
        remove_element(alleles, fasta[i])
        if n>1:
            variable.update({i+1: alleles})
            if verbose:
                print(i+1, j+1, fasta[i], alleles, sep='\t')
    return variable


def get_transcript_variant_number(transcript_record):
    '''
    Get transcript variant number from sequence definition.

    Return 0 if there is no indication.

    Parameters
    ----------
    transcript_record : dict
        Transcript record.

    Return
    ------
    int
        Transcript variant number.
    '''
    gene_def = transcript_record['GBSeq_definition']
    gene_def_list1 = gene_def.split('transcript variant ')
    if len(gene_def_list1) == 1:
        return 0
    else:
        gene_def_list2 = gene_def_list1[1].split(', ')
        return int(gene_def_list2[0])


def get_gene_abbrev(transcript_record):
    '''
    Get abbreviated gene name from sequence definition.

    Parameters
    ----------
    transcript_record : dict
        Transcript record.

    Return
    ------
    str
        Abbreviated gene name.
    '''
    gene_def = transcript_record['GBSeq_definition']
    abbrev1 = gene_def.split('(')
    abbrev2 = abbrev1[1].split(')')
    return abbrev2[0]


def get_CDS(transcript_record):
    '''
    Get the transcript, CDS, and protein sequences from a transcript record.

    Parameters
    ----------
    transcript_record : dict
        Transcript record.

    Return
    ------
    tuple (beg, end, seq, cds_seq, prot_seq)
        beg : int
            Beginning of CDS in transcript sequence (from 0).
        end : int
            End of CDS in transcript sequence (from 0).
        seq : Seq object
            Transcript sequence.
        cds_seq : Seq object
            Coding sequence.
        prot_seq : Seq object
            Protein sequence.
    '''
    seq = transcript_record['GBSeq_sequence']
    nfeatures = len(transcript_record['GBSeq_feature-table'])
    for i in range(nfeatures):
        f = transcript_record['GBSeq_feature-table'][i]
        if f['GBFeature_key'] == 'CDS':
            loc = f['GBFeature_location'].split('..')
            beg = int(loc[0]) - 1
            end = int(loc[1])
            nquals = len(f['GBFeature_quals'])
            for j in range(nquals):
                q = f['GBFeature_quals'][j]
                if q['GBQualifier_name'] == 'translation':
                    prot_seq = Seq(q['GBQualifier_value'])
                    cds_seq = Seq(seq[beg:end])
                    return beg, end, seq, cds_seq, prot_seq


def mutate(seq, i, new):
    '''
    Mutate DNA sequence at a particular site to a particular nucleotide.

    Parameters
    ----------
    seq : str or Seq object.
        DNA sequence.
    i : int
        Site number (from 0).
    new : str
        Mutant nucleotide.

    Returns
    -------
    str
        Mutant sequence.
    '''
    L = len(seq)
    assert i < L
    assert seq[i] != new
    return seq[:i] + new + seq[(i+1):]


def get_all_site_mutants(seq, i):
    '''
    Mutate DNA sequence at a particular site to all alternative nucleotides.

    Parameters
    ----------
    seq : str or Seq object.
        DNA sequence.
    i : int
        Site number (from 0).

    Returns
    -------
    list
        Mutant sequences.
    '''
    mutants = []
    for j in nucl:
        if seq[i] != j:
            mutants.append(mutate(seq, i, j))
    return mutants


def get_all_mutants(seq):
    '''
    Mutate DNA sequence at each site to all alternative nucleotides.

    Parameters
    ----------
    seq : str or Seq object.
        DNA sequence.

    Returns
    -------
    list
        Mutant sequences.
    '''
    L = len(seq)
    mutants = []
    for i in range(L):
        mutants.extend(get_all_site_mutants(seq, i))
    return mutants


def get_codon_mutations(codon, verbose=False):
    '''
    Generate all single mutants of a codon and evaluate their phenotypic effects.

    Parameters
    ----------
    codon : Seq object
        Codon DNA sequence.
    verbose : bool
        Whether to show breakdown of individual sequences.

    Returns
    -------
    tuple of ints
        Counts of synonymous, nonsynonymous, and nonsense mutations.
    '''
    assert len(codon) == 3
    aa = codon.translate()
    if verbose:
        print('Reference:  ', codon, aa)
    syn = 0
    nonsyn = 0
    nonsen = 0
    mut_codons = get_all_mutants(codon)
    for mut_codon in mut_codons:
        mut_aa = mut_codon.translate()
        if mut_aa == aa:
            syn += 1
        else:
            if mut_aa == '*':
                nonsen += 1
            else:
                nonsyn += 1
        if verbose:
            print('   Mutant:  ', mut_codon, mut_aa, syn, nonsyn, nonsen)
    return syn, nonsyn, nonsen


def get_mutation_counts(seq, verbose=False):
    '''
    Generate all single mutants of a sequence and evaluate their phenotypic effects.

    Parameters
    ----------
    seq : Seq object
        DNA sequence.
    verbose : bool
        Whether to show breakdown of individual codons.

    Returns
    -------
    tuple of ints
        Counts of synonymous, nonsynonymous, and nonsense mutations.
    '''
    syn = 0
    nonsyn = 0
    nonsen = 0
    L = len(seq.translate(cds=True))
    for i in range(L):
        codon = seq[(i * 3):((i + 1) * 3)]
        newsyn, newnonsyn, newnonsen = get_codon_mutations(codon, verbose)
        syn += newsyn
        nonsyn += newnonsyn
        nonsen += newnonsen
    return syn, nonsyn, nonsen


########################
### ClinVar variants ###
########################


def get_variant_ids(gene):
    '''
    Get variant IDs from gene name.

    Parameters
    ----------
    gene : str
        Gene name.

    Return
    ------
    list
        ClinVar IDs.

    '''
    search_term = "{0}[Gene Name]".format(gene)
    handle = Entrez.esearch(db="clinvar", term=search_term, retmax=10000)
    search = Entrez.read(handle)
    search_ids = search['IdList']
    handle.close()
    return search_ids


def get_variant(variant):
    '''
    Get ClinVar record from ID.

    Parameters
    ----------
    variant : str
        ClinVar record ID.

    Return
    ------
    dict
        ClinVar record.
    '''
    handle = Entrez.efetch(db="clinvar", id=variant, rettype='vcv', is_varationid=True, from_esearch=True)
    tmp = x2d.parse(handle.read().decode('utf-8'))
    record = tmp['ClinVarResult-Set']['VariationArchive']
    handle.close()
    return record


def get_variant_details(variant_record):
    '''
    Extract information from the @VariationName field.

    Parameters
    ----------
    variant_record : dict
        ClinVar record.

    Return
    ------
    tuple (gene_name, transcript, nucl, prot, mol_conseq)
        gene_name : str
            Abbreviated gene name.
        transcript : str
            Transcript name.
        nucl : str
            Nucleotide change.
        prot : str
            Protein change.
        mol_conseq : str
            Molecular consequence.
        effect : str
            Phenotypic effect.
    '''
    var_name = variant_record['@VariationName']
    var_name_split1 = var_name.split('(')
    var_name_split2 = var_name_split1[1].split(')')
    gene_name = var_name_split2[0]
    transcript = var_name_split1[0]
    nucl = var_name_split2[1][3:-1]
    if len(var_name_split1) == 3:
        prot = var_name_split1[2][2:-1]
        if prot[-1] == '=':
            mol_conseq = 'synonymous'
        elif prot[-3:] == 'Ter':
            mol_conseq = 'nonsense'
        else:
            mol_conseq = 'missense'
    else:
        prot = 'none'
        mol_conseq = 'other'
    if variant_record['@RecordType'] == 'interpreted':
        rec_type = 'InterpretedRecord'
    elif variant_record['@RecordType'] == 'included':
        rec_type = 'IncludedRecord'
    effect = variant_record[rec_type]['Interpretations']['Interpretation']['Description']
    return gene_name, transcript, nucl, prot, mol_conseq, effect


def get_species(record):
    '''
    Retrieve species name from NCBI sequence record.

    Parameters
    ----------
    record : Bio.SeqRecord
        Sequence record.

    Returns
    -------
    str
        Species name.
    '''
    description = record.description
    species = description.split(' [')[1][:-1]
    return species


def get_alleles(gene, protein, site, verbose=True):
    """
    Get list of alleles at a particular site in the alignment.

    Parameters
    ----------
    gene : str
        Gene name.
    protein : str
        NCBI record ID.
    site : int
        Amino acid site.
    verbose : type
        Description of parameter `verbose`.
    """
    for record in SeqIO.parse("fasta/{0}.fasta".format(gene), "fasta"):
        recordid = record.description.split(' ')[0]
        if recordid == protein:
            fasta = record
            break
    align = AlignIO.read('aln/{0}.aln'.format(gene), 'fasta')
    for record in align:
        if record.id == protein:
            aln = record
            break
    i = site - 1
    alni = get_aln_positions(gene, protein)
    j = alni[i]
    if verbose:
        print('     Gene:', gene)
        print('    Human:', i+1, fasta.seq[i])
        print('Alignment:', j+1, aln.seq[j])
    counts = pd.Series(list(align[:, j])).value_counts()
    if len(counts) == 1:
        return False, i+1, fasta.seq[i], counts
    else:
        if verbose:
            print(counts)
            print()
            for record in align:
                print(record[j], record.id, get_species(record))
        return True, i+1, fasta.seq[i], counts


def get_species_with_allele(gene, protein, site, allele):
    """
    Get species with a particular amino acid at a particular site in the
    alignment.

    Parameters
    ----------
    gene : str
        Gene name.
    protein : str
        NCBI record ID.
    site : int
        Amino acid site.
    allele : int
        Particular amino acid at site.

    Returns
    -------
    list
        (species, protein ID) : tup
    """
    for record in SeqIO.parse("fasta/{0}.fasta".format(gene), "fasta"):
        recordid = record.description.split(' ')[0]
        if recordid == protein:
            fasta = record
            break
    align = AlignIO.read('aln/{0}.aln'.format(gene), 'fasta')
    for record in align:
        if record.id == protein:
            aln = record
            break
    i = site - 1
    alni = get_aln_positions(gene, protein)
    j = alni[i]
    alleles = pd.Series(list(align[:, j]))
    ii = alleles[alleles==allele].index.tolist()
    spp = []
    for i in ii:
        spp.append((get_species(align[i]), align[i].id))
    return spp


def local_compare_to_human(gene, human, site, nonhuman):
    for record in SeqIO.parse("fasta/{0}.fasta".format(gene), "fasta"):
        recordid = record.description.split(' ')[0]
        if recordid == human:
            fasta = record
            break
    align = AlignIO.read('aln/{0}.aln'.format(gene), 'fasta')
    L = len(align[0].seq)
    for record in align:
        if record.id == human:
            aln = record
        elif record.id == nonhuman:
            nonaln = record
    i = site - 1
    alni = get_aln_positions(gene, human)
    j = alni[i]
    diff = 0
    out = 0
    gaps = 0
    for delta in range(-10, 11):
        if (0 <= j+delta < L):
            print(delta, aln.seq[j+delta], nonaln.seq[j+delta])
            if (delta != 0) and (aln.seq[j+delta] != nonaln.seq[j+delta]):
                diff += 1
                if (aln.seq[j+delta] == '-') or (nonaln.seq[j+delta] == '-'):
                    gaps += 1
        else:
            print(delta, 'outside sequence')
            out += 1
    return diff, gaps, out


def global_compare_to_human(gene, human, nonhuman, verbose):
    differences = {}
    for record in SeqIO.parse("fasta/{0}.fasta".format(gene), "fasta"):
        recordid = record.description.split(' ')[0]
        if recordid == human:
            fasta = record
            break
    align = AlignIO.read('aln/{0}.aln'.format(gene), 'fasta')
    L = len(align[0].seq)
    for record in align:
        if record.id == human:
            aln = record
        elif record.id == nonhuman:
            nonaln = record
    for i in range(L):
        if (aln.seq[i] != nonaln.seq[i]) and (aln.seq[i] != '-') and (nonaln.seq[i] != '-'):
            if verbose:
                print(i+1, aln.seq[i], nonaln.seq[i])
            differences.update({i+1: (aln.seq[i], nonaln.seq[i])})
    return differences
