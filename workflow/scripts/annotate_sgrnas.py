import logging
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import Align
import pyranges as pr
import pyfaidx as px
from tqdm import tqdm
import concurrent.futures

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

def process_sgrna(sgrna):
    gene = sgrna.split("_")[0]
    gene_gtf = GTF[GTF["gene_name"] == gene]
    if len(gene_gtf) == 0:
        return (gene, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, True)
    
    chrom = gene_gtf["Chromosome"].values[0]
    strand = gene_gtf["Strand"].values[0]
    seq = sgrnas[sgrna][0:20].seq
    
    gene_start = gene_gtf[gene_gtf["Feature"] == "gene"]["Start"].values
    assert len(gene_start) == 1, f"Gene {gene} has multiple start positions: {gene_start}"
    gene_start = gene_start[0]
    gene_end = gene_gtf[gene_gtf["Feature"] == "gene"]["End"].values
    assert len(gene_end) == 1, f"Gene {gene} has multiple end positions: {gene_end}"
    gene_end = gene_end[0]
    gene_seq = FASTA[chrom][gene_start:gene_end].seq
    
    sg_start = gene_seq.find(seq)
    if (sg_start == -1):
        gene_seq = str(Seq(gene_seq).reverse_complement())
        sg_start = gene_seq.find(seq)
        if sg_start == -1:
            raise ValueError(f"Could not find sgRNA sequence {seq} in gene {gene}")
    sg_end = sg_start + len(seq)
    sg_start = gene_start + sg_start + 1
    sg_end = gene_start + sg_end + 1
    
    #return (gene, chrom, sg_start, sg_end, strand, np.nan, seq)
    return (gene, chrom, sg_start, sg_end, strand, np.nan, seq, False)

# Set up parallel processing
threads = snakemake.threads
with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
    sgrna_keys = list(sgrnas.keys())
    results = list(tqdm(executor.map(process_sgrna, sgrna_keys), total=len(sgrna_keys)))

for result in results:
    gene, chrom, sg_start, sg_end, strand, exon_number, seq = result
    gene, chrom, sg_start, sg_end, strand, exon_number, seq, gene_not_found = result
    if gene_not_found:
        genes_not_found.append(gene)
    chromosomes.append(chrom)
    sg_starts.append(sg_start)
    sg_ends.append(sg_end)
    strands.append(strand)
    exon_numbers.append(exon_number)
    sequences.append(seq)
    
           
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
    