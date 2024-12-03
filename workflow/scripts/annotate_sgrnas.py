import logging
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import pyranges as pr
import pyfaidx as px
from tqdm import tqdm

"""
Annotate each sgRNA with:
- gene
- chromosome
- start
- end
- strand
- exon number
- sequence
"""

log_file = snakemake.log[0]
logging.basicConfig(format='%(levelname)s:%(message)s', 
                    level=logging.DEBUG,
                    handlers=[logging.FileHandler(log_file)])

sgrnas = px.Fasta(snakemake.params["sgrnas"])
GTF = pr.read_gtf(snakemake.input["gtf"]).df
FASTA = px.Fasta(snakemake.input["fasta"])

# Lists to store annotations
genes = []
chromosomes = []
sg_starts = []
sg_ends = []
strands = []
exon_numbers = []
sequences = []

genes_not_found = []

for sgrna in tqdm(sgrnas.keys()):
    # Get gene name
    gene = sgrna.split("_")[0]
    genes.append(gene)
    
    # Subset GTF for gene
    gene_gtf = GTF[GTF["gene_name"] == gene]
    if len(gene_gtf) == 0:
        genes_not_found.append(gene)
        # Add NaN values to all lists
        chromosomes.append(np.nan)
        strands.append(np.nan)
        sequences.append(np.nan)
        sg_starts.append(np.nan)
        sg_ends.append(np.nan)
        exon_numbers.append(np.nan)
        continue
        
    #assert len(gene_gtf) > 0, f"Gene {gene} not found in GTF"
    
    # Get chromosome
    chrom = gene_gtf["Chromosome"].values[0]
    chromosomes.append(chrom)
    
    # Get strand
    strand = gene_gtf["Strand"].values[0]
    strands.append(strand)
    
    # Get sgRNA sequence
    seq = sgrnas[sgrna][0:20].seq
    sequences.append(seq)
        
    # Get start and end of sgRNA sequence inside gene
    gene_start = gene_gtf[gene_gtf["Feature"] == "gene"]["Start"].values
    assert len(gene_start) == 1, f"Gene {gene} has multiple start positions: {gene_start}"
    gene_start = gene_start[0]
    gene_end = gene_gtf[gene_gtf["Feature"] == "gene"]["End"].values
    assert len(gene_end) == 1, f"Gene {gene} has multiple end positions: {gene_end}"
    gene_end = gene_end[0]
    gene_seq = FASTA[chrom][gene_start:gene_end].seq
    
    # Find start and end of sgRNA sequence inside gene
    sg_start = gene_seq.find(seq)
        
    if sg_start == -1:
        #logging.warning(f"Could not find sgRNA sequence {seq} in gene {gene}")
        seq = str(Seq(seq).reverse_complement())
        sg_start = gene_seq.find(seq)
        if sg_start == -1:
            raise ValueError(f"Could not find sgRNA sequence {seq} in gene {gene}")
    sg_end = sg_start + len(seq)
    # Convert coordinates inside gene to genome coordinates
    # Also convert to 1-based coordinates
    sg_start = gene_start + sg_start + 1
    sg_starts.append(sg_start)
    sg_end = gene_start + sg_end + 1
    sg_ends.append(sg_end)
    
    # Get exon number sgRNA is in:
    # Get all exons of gene
    # Check which exon overlaps with sgRNA targeting coordinates
    # NOTE: sgRNA may also be part in intron
    exons = GTF[GTF["gene_name"] == gene]
    exons = exons[exons["Feature"] == "exon"]
    
    # Check for sgRNA-exon overlap within exon 
    exons_within = exons[(exons["Start"] <= sg_start) & (exons["End"] >= sg_end)]
    
    # Check for sgRNA-exon partial overlap upstream
    exons_partial_upstream = exons[(exons["Start"] >= sg_start) & ((exons["End"] >= sg_end) & (exons["Start"] <= sg_end))] 
    
    # Check for sgRNA-exon partial overlap downstream
    exons_partial_downstream = exons[(exons["Start"] <= sg_start) & ((exons["End"] <= sg_end) & (exons["End"] >= sg_start))]
    
    # Combine all exons
    all_exons = pd.concat([exons_within, exons_partial_upstream, exons_partial_downstream])
    found_exons = all_exons["exon_number"].values
    if len(all_exons["exon_id"].values) == 1:
        exon_numbers.append(found_exons[0])
    else:
        # Find canonical exon
        canonical_exon = all_exons[all_exons["tag"] == "Ensembl_canonical"]
        if len(canonical_exon) == 1:
            exon_numbers.append(canonical_exon["exon_number"].values[0])
        else:
            # Check if all exon_numbers are identical
            if len(set(found_exons)) == 1:
                exon_numbers.append(found_exons[0])
            else:
                # Just pick the first exon
                #logging.warning(f"sgRNA {sgrna} targets multiple exons: {found_exons} (picked first)")
                #exon_numbers.append(found_exons[0])
                raise ValueError(f"sgRNA {sgrna} targets multiple exons: {found_exons}")
       
# Create data frame
df = pd.DataFrame({
    "CODE": list(sgrnas.keys()),
    "GENES": genes,
    "CHRM": chromosomes,
    "STARTpos": sg_starts,
    "ENDpos": sg_ends,
    "STRAND": strands,
    "EXONE": exon_numbers,
    "seq": sequences
})   
    
# Save data frame
df.to_csv(snakemake.output[0], sep="\t", index=False)

# Log genes not found
if len(genes_not_found) > 0:
    genes_not_found = "\n".join(genes_not_found)
    logging.warning(f"Genes not found in GTF (most likely due to wrong synonym used...): {genes_not_found}")   
    