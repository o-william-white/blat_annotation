import argparse
import subprocess
from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import os

# argparse
parser = argparse.ArgumentParser()
parser.add_argument("--query",      help = "Input fasta to annotate",        required=True)
parser.add_argument("--reference",  help = "NCBI accession for reference",   required=True)
parser.add_argument("--email",      help = "Email for Entrez search",        required=True)
parser.add_argument("--prefix",     help = "Output prefix",                  required=False, default="blat_annotation")
parser.add_argument("--circular",   help = "Input query is circular",        required=False, action="store_true")
parser.add_argument("--keep",       help = "Keep reference and blat files",  required=False, action="store_true")
args = parser.parse_args()

# define user email 
Entrez.email = args.email

# read in query
query_record = list(SeqIO.parse(args.query, "fasta"))[0]

# if the query is circular, append the first 10Kb to the end of the sequence query
if args.circular:
    print("Treating query as circular")
    query_record.seq = query_record.seq + query_record.seq[0:10000]

# download genbank of reference data and write annotated genes to fasta

# efetch reference accession
handle = Entrez.efetch(db="nucleotide", id=args.reference, rettype="gb", retmode="text")
reference_record = SeqIO.read(handle, "gb")
handle.close()

# list and dictionary to store annotations and counts respectively
list_coding, list_non_coding, dict_genes, dict_counts = [], [], {}, {}

# loop through sequence feature and append to list
for feature in reference_record.features:
    # cds features
    if (feature.type in ["CDS", "tRNA", "rRNA"]) and ("gene" in feature.qualifiers.keys()):
        # get annotation info
        gene_name = str(feature.qualifiers["gene"][0])
        gene_seq = feature.extract(reference_record)
        gene_seq.description = ""
        if dict_genes.get(gene_name) == None:
            dict_genes[gene_name] = str(gene_seq.seq)
            dict_counts[gene_name] = 1
        elif dict_genes[gene_name] == str(gene_seq.seq): 
            print(f"Ignoring duplicate annotation for {gene_name}")
            continue
        elif dict_genes[gene_name] != str(gene_seq.seq):
            dict_counts[gene_name] += 1
            print(f"Note {dict_counts[gene_name]} sequences have been extracted under the name {gene_name}")
        if feature.type == "CDS": 
            gene_seq.id = f"cds-blatx_{gene_name}_{dict_counts[gene_name]};gbkey={feature.type};gene={gene_name}"
            list_coding.append(gene_seq)
        elif (feature.type == "tRNA") or (feature.type == "rRNA"):
            gene_seq.id = f"cds-blatn_{gene_name}_{dict_counts[gene_name]};gbkey={feature.type};gene={gene_name}"
            list_non_coding.append(gene_seq)

# write reference fasta
with open(f"{args.prefix}_reference_coding.fasta", "w") as coding_fas:
    SeqIO.write(list_coding, coding_fas, "fasta")

with open(f"{args.prefix}_reference_non_coding.fasta", "w") as non_coding_fas:
    SeqIO.write(list_non_coding, non_coding_fas, "fasta")

# run blat
subprocess.run(["blat", f"{args.prefix}_reference_coding.fasta",      args.query, "-t=dnax", "-q=dnax", "-out=blast8", f"{args.prefix}_output_coding.txt"])
subprocess.run(["blat", f"{args.prefix}_reference_non_coding.fasta",  args.query, "-t=dna",  "-q=dna",  "-out=blast8", f"{args.prefix}_output_non_coding.txt"])

# read blat output into pandas dataframes
blast8_columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
blat_coding     = pd.read_csv(f"{args.prefix}_output_coding.txt",     sep="\t", names=blast8_columns)
blat_non_coding = pd.read_csv(f"{args.prefix}_output_non_coding.txt", sep="\t", names=blast8_columns)

# select top hits based on bitscores
blat_coding_top     = blat_coding.loc[blat_coding.groupby('sseqid')['bitscore'].idxmax()]
blat_non_coding_top = blat_non_coding.loc[blat_non_coding.groupby('sseqid')['bitscore'].idxmax()]

# concat coding and noncoding blat outputs
dat = pd.concat([blat_coding_top, blat_non_coding_top], axis=0)
dat = dat.sort_values(by="qstart").reset_index(drop=True)

# create bed
bed = dat[["qseqid", "qstart", "qend", "sseqid"]].rename(columns={
    "qseqid" : "chrom",
    "qstart" : "start",
    "qend" : "end",
    "sseqid" : "name"
})

# write bed
bed.to_csv(f"{args.prefix}_example.bed", sep="\t", header=True, index=False)

# remove intermediate files
if not args.keep:
    os.remove(f"{args.prefix}_output_coding.txt")
    os.remove(f"{args.prefix}_output_non_coding.txt")
    os.remove(f"{args.prefix}_reference_coding.fasta")
    os.remove(f"{args.prefix}_reference_non_coding.fasta")
