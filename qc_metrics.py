import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from Bio import SeqIO

def count_genes(bam_file):
    gene_counts = {}
    print(bam_file)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        lengths = []
        for read in bam:
            lengths.append(read.query_length)
            gene_name = read.reference_name  # Using reference_name as a proxy for gene name
            if gene_name not in gene_counts:
                gene_counts[gene_name] = 0
            gene_counts[gene_name] += 1
    return gene_counts , lengths

def exploring_bam(location):
    df = pd.read_csv(location,sep='\t',header=0)
    df[['cell_line', 'condition', 'replicate']] = df['sample_id'].str.split('_', expand=True)
    multi_index_df = df.set_index(['cell_line', 'condition', 'replicate'])

    gene_counts_df = pd.DataFrame()
    read_lengths_df = pd.DataFrame()

    for index, row in multi_index_df.iterrows():
        counts, lengths = count_genes(row['bam'])
        counts_df = pd.DataFrame(counts, index=[index])
        counts_df["lengths"] = [lengths]
        gene_counts_df = pd.concat([gene_counts_df, counts_df], axis=0)

    gene_counts_df = gene_counts_df.fillna(0)
    columns_to_sum = gene_counts_df.columns.difference(['lengths'])
    gene_counts_df_percentage = gene_counts_df[columns_to_sum].div(gene_counts_df[columns_to_sum].sum(axis=1), axis=0) * 100

    ax = gene_counts_df.plot(kind='bar', stacked=True, figsize=(12, 8), legend=False)
    ax.set_xlabel('Genes')
    ax.set_ylabel('Counts')
    ax.set_title('')
    plt.tight_layout()
    plt.savefig('gene_counts_comparison_counts.png')
    plt.close()

    fig, ax = plt.subplots(figsize=(12, 8))
    lengths_df = gene_counts_df["lengths"]
    lengths_flat = lengths_df.iloc[1]

    length_counts = pd.Series(lengths_flat).value_counts().sort_index()
    plt.plot(length_counts.index, length_counts.values, marker='o', linestyle='-', color='b')
    plt.xlabel('Lengths')
    plt.ylabel('Counts')
    plt.title('')
    plt.savefig("read_length.png")
    plt.close()

def count_unique_mapped_reads(bam_file_path):
    input_bam = pysam.AlignmentFile(bam_file_path, "rb")
    unique_mapped_count = 0
    for read in input_bam:
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary and read.mapping_quality == 60:
            unique_mapped_count += 1
    input_bam.close()
    return unique_mapped_count

def calculate_average_phred_score(qualities):
    return sum(qualities) / len(qualities)

def process_fastq_file(file_path):
    """Read a FASTQ file and return a list of tuples (read_length, average_phred_score)."""
    results = []

    with open(file_path, 'r') as fastq_file:
        for record in SeqIO.parse(fastq_file, 'fastq'):
            read_length = len(record.seq)
            avg_phred_score = calculate_average_phred_score(record.letter_annotations["phred_quality"])
            results.append((read_length, avg_phred_score))
    return results

def create_scatterplot(scatterplot_list,data_type,threshold):
    """Create a scatterplot from the list of tuples."""
    read_lengths, avg_phred_scores = zip(*scatterplot_list)
    avg_phred_scores = np.array(avg_phred_scores)  # Convert to NumPy array

    plt.figure(figsize=(10, 6))
    
    # Scatter plot
    plt.scatter(read_lengths, avg_phred_scores, alpha=0.5, label='Data Points')
    # Shaded red area from 0 to y-coordinate 20
    plt.xlim(0,60000)
    plt.ylim(0,100)
    plt.title('')
    plt.xlabel('Read Length')
    plt.ylabel('Average Phred Score')
    plt.axhline(y=threshold, color='red', linestyle='--', label=f"Threshold (y={threshold})")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{data_type}_scatterplot.png", format="png")
    plt.close()

def create_histogram(scatterplot_list,data_type):
    bin_ranges = [0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40]
    read_lengths, avg_phred_scores = zip(*scatterplot_list)
    rs = pd.Series(read_lengths)
    ps = pd.Series(avg_phred_scores)
    print(data_type,rs.describe(),ps.describe())
    plt.figure(figsize=(10, 6))
    plt.hist(avg_phred_scores, bins=bin_ranges, edgecolor='black')
    plt.ylim(0,2000000)
    plt.title('')
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.savefig(f"{data_type}_Phred_score_histogram.png", format="png")

def read_analysis(location,threshold):
    df = pd.read_csv(location,sep='\t',header=0)
    # Create an empty column to store the results in the DataFrame
    df["read_lengths"] = None
    df["avg_phred_scores"] = None

    # Process each row in the DataFrame
    scatterplot_list = []
    for index, row in df.iterrows():
        results = process_fastq_file(row["processed_fastq"])
        read_lengths, avg_phred_scores = zip(*results)
        df.at[index,"processed_reads"] = len(read_lengths)
        df.at[index,"mapped_reads"] = count_unique_mapped_reads(row["bam"])
        
        scatterplot_list.extend(results)

    create_scatterplot(scatterplot_list,"processed",threshold)
    create_histogram(scatterplot_list,"processed")

    scatterplot_list = []
    for index, row in df.iterrows():
        results = process_fastq_file(row["raw_fastq"])
        read_lengths, avg_phred_scores = zip(*results)
        df.at[index,"basecalled_reads"] = len(read_lengths)
        df.at[index, "read_lengths"] = read_lengths
        df.at[index, "avg_phred_scores"] = avg_phred_scores
        scatterplot_list.extend(results)

    df.to_csv("sequenced_amounts.tsv",sep='\t',columns=["sample_id","basecalled_reads","processed_reads","mapped_reads"])

    create_scatterplot(scatterplot_list,"raw",threshold)
    create_histogram(scatterplot_list,"raw")

    plt.figure(figsize=(10, 6))
    sns.boxplot(x="sample_id", y="read_lengths", data=df.explode("read_lengths"), showfliers=False)
    plt.title('')
    plt.xlabel('Sample ID')
    plt.ylabel('Read Length')
    plt.xticks(rotation=90)
    plt.tight_layout() 
    plt.savefig("raw_boxplot_read_lengths_no_outliers.png", format="png")
    plt.close()

    # Create a boxplot for average Phred scores with removed outliers and rotated labels
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="sample_id", y="avg_phred_scores", data=df.explode("avg_phred_scores"), showfliers=False)
    plt.title('')
    plt.xlabel('Sample ID')
    plt.ylabel('Average Phred Score')
    plt.xticks(rotation=90)
    plt.tight_layout() 
    plt.savefig("raw_boxplot_avg_phred_scores_no_outliers.png", format="png")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read and display a TSV file.')
    parser.add_argument('file_path', metavar='path', type=str, help='Path to the TSV file')
    parser.add_argument('threshold', type=int, help='value of phred cutoff')

    args = parser.parse_args()

    # Call the function to read and display the TSV file
    exploring_bam(args.file_path)
    read_analysis(args.file_path,args.threshold)
    print("QC metrics complete.")
