import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import seaborn as sns
import numpy as np
import pysam

# Load the TSV file containing sample metadata and file paths

def count_unique_mapped_reads(bam_file_path):
    """
    Count the number of uniquely mapped reads in a BAM file.

    Args:
        bam_file_path (str): Path to the BAM file.

    Returns:
        int: Count of uniquely mapped reads.
    """
    input_bam = pysam.AlignmentFile(bam_file_path, "rb")
    unique_mapped_count = 0
    for read in input_bam:
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary and read.mapping_quality == 60:
            unique_mapped_count += 1
    input_bam.close()
    return unique_mapped_count

def calculate_average_phred_score(qualities):
    """
    Calculate the average Phred quality score for a list of quality scores.

    Args:
        qualities (list): List of Phred quality scores.

    Returns:
        float: Average Phred quality score.
    """
    return sum(qualities) / len(qualities)

def process_fastq_file(file_path):
    """
    Process a FASTQ file to extract read lengths and average Phred scores.

    Args:
        file_path (str): Path to the FASTQ file.

    Returns:
        list: List of tuples containing read lengths and average Phred scores.
    """
    results = []

    with open(file_path, 'r') as fastq_file:
        for record in SeqIO.parse(fastq_file, 'fastq'):
            read_length = len(record.seq)
            avg_phred_score = calculate_average_phred_score(record.letter_annotations["phred_quality"])
            results.append((read_length, avg_phred_score))
    return results

def create_scatterplot(scatterplot_list, data_type):
    """
    Create a scatterplot of read lengths vs. average Phred scores.

    Args:
        scatterplot_list (list): List of tuples (read_length, average_phred_score).
        data_type (str): Type of data being plotted (e.g., 'processed', 'raw').

    Returns:
        None
    """
    read_lengths, avg_phred_scores = zip(*scatterplot_list)
    avg_phred_scores = np.array(avg_phred_scores)  # Convert to NumPy array

    plt.figure(figsize=(10, 6))
    
    # Scatter plot
    plt.scatter(read_lengths, avg_phred_scores, alpha=0.5, label='Data Points')
    # Shaded red area from 0 to y-coordinate 20
    plt.xlim(0,60000)
    plt.ylim(0,100)
    plt.title('Read Length vs Average Phred Score')
    plt.xlabel('Read Length')
    plt.ylabel('Average Phred Score')
    plt.axhline(y=16, color='red', linestyle='--', label='Threshold (y=16)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{data_type}_scatterplot.png", format="png")
    plt.close()

def create_histogram(scatterplot_list, data_type):
    """
    Create a histogram of average Phred scores.

    Args:
        scatterplot_list (list): List of tuples (read_length, average_phred_score).
        data_type (str): Type of data being plotted (e.g., 'processed', 'raw').

    Returns:
        None
    """
    bin_ranges = [0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40]
    read_lengths, avg_phred_scores = zip(*scatterplot_list)
    plt.figure(figsize=(10, 6))
    plt.hist(avg_phred_scores, bins=bin_ranges, edgecolor='black')
    plt.ylim(0,2000000)
    plt.title('Histogram')
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.savefig(f"{data_type}_Phred_score_histogram.png", format="png")

# Add columns to the DataFrame to store read lengths and average Phred scores
df["read_lengths"] = None
df["avg_phred_scores"] = None

# Process each row in the DataFrame to analyze processed FASTQ files
scatterplot_list = []
for index, row in df.iterrows():
    # Extract read lengths and average Phred scores from processed FASTQ files
    results = process_fastq_file(row["processed_fastq"])
    read_lengths, avg_phred_scores = zip(*results)
    df.at[index, "processed_reads"] = len(read_lengths)
    df.at[index, "mapped_reads"] = count_unique_mapped_reads(row["bam"])
    scatterplot_list.extend(results)

# Generate scatterplot and histogram for processed data
create_scatterplot(scatterplot_list, "processed")
create_histogram(scatterplot_list, "processed")

# Process each row in the DataFrame to analyze raw FASTQ files
scatterplot_list = []
for index, row in df.iterrows():
    # Extract read lengths and average Phred scores from raw FASTQ files
    results = process_fastq_file(row["fastq"])
    read_lengths, avg_phred_scores = zip(*results)
    df.at[index, "basecalled_reads"] = len(read_lengths)
    df.at[index, "read_lengths"] = read_lengths
    df.at[index, "avg_phred_scores"] = avg_phred_scores
    scatterplot_list.extend(results)

# Save the results to a TSV file
df.to_csv("sequenced_amounts.tsv", sep='\t', columns=["alias", "basecalled_reads", "processed_reads", "mapped_reads"])

# Generate scatterplot and histogram for raw data
create_scatterplot(scatterplot_list, "raw")
create_histogram(scatterplot_list, "raw")

# Create a boxplot for read lengths with removed outliers
plt.figure(figsize=(10, 6))
sns.boxplot(x="alias", y="read_lengths", data=df.explode("read_lengths"), showfliers=False)
plt.title('Boxplot of Read Lengths for Each Sample')
plt.xlabel('Sample ID')
plt.ylabel('Read Length')
plt.xticks(rotation=90)
plt.grid(True)
plt.tight_layout()
plt.savefig("raw_boxplot_read_lengths_no_outliers.png", format="png")
plt.show()

# Create a boxplot for average Phred scores with removed outliers
plt.figure(figsize=(10, 6))
sns.boxplot(x="alias", y="avg_phred_scores", data=df.explode("avg_phred_scores"), showfliers=False)
plt.title('Boxplot of Average Phred Scores for Each Sample')
plt.xlabel('Sample ID')
plt.ylabel('Average Phred Score')
plt.xticks(rotation=90)
plt.grid(True)
plt.tight_layout()
plt.savefig("raw_boxplot_avg_phred_scores_no_outliers.png", format="png")
plt.show()

