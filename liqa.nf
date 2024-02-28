#!/usr/bin/env nextflow

params.fast5 = "{pathway}/pod5"
params.sample_sheet = "{pathway}/liqa_reference.tsv"
params.phred_score = 10
params.q_score = 60
params.outdir = "liqa.output"
params.basecalling = false
params.fastq_dir = "{pathway}/guppy.combined.out/pass"
params.reference_gtf = "{pathway}/Homo_sapiens_GRCh38_110.gtf"
params.reference_fna = "{pathway}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.reference_bed = "{pathway}/Homo_sapiens_GRCh38_110.bed"

process basecalling {
    label 'basecalling'

    publishDir "${params.outdir}/basecalling", mode: 'copy'

    input:
        path fast5

    output:
        path "guppy.out/pass"
    
    script:
    """
    guppy_basecaller \
	--barcode_kits "SQK-NBD114-24" \
	--disable_pings -i ${fast5} -s guppy.out \
	-c "dna_r10.4.1_e8.2_400bps_hac.cfg" \
	--recursive -x "auto" \
	--num_callers ${task.cpus} \
	--gpu_runners_per_device 4 \
	--chunks_per_runner 1024 \
	--chunk_size 1000 \
	--chunks_per_caller 10000 \
	--compress_fastq
    """
}


process makeRefgene{
    executor="local"

    input:
        path reference_gtf

    output:
        path "reference.refgene"

    script:
    """
    liqa -task refgene \
    -ref ${reference_gtf} \
    -format gtf \
    -out reference.refgene
    """
}

process removeAdaptors {
    label "remove_adaptor"

    publishDir "${params.outdir}/raw_reads", mode: 'copy', pattern: "${sample_id}.combined.fastq"

    input:
        tuple val(sample_id) , val(barcode) , val(cell_line), val(condition), val(replicate), path(fastq_dir)

    output:
        tuple val(sample_id) , val(barcode) , val(cell_line), val(condition), val(replicate), path("trimmed.fastq"), path("locations.tsv"), emit: reads
        path "${sample_id}.combined.fastq"

    script:
    """
    if [ -d "${fastq_dir}/pass" ]; then
        fastq_dir="${fastq_dir}/pass"
    else
        fastq_dir="${fastq_dir}"
    fi

    cat \$fastq_dir/*.fastq.gz | gunzip -c > ${sample_id}.combined.fastq

    porechop -i "${sample_id}.combined.fastq" -o "trimmed.fastq" --threads ${task.cpus}


    echo -ne "${sample_id}\t${cell_line}\t${condition}\t${replicate}\t${baseDir}/${params.outdir}/raw_reads/${sample_id}.combined.fastq" > "locations.tsv"
    """
}

process preprocessingReads {
    label "preprocessing"

    publishDir "${params.outdir}/preprocessed", mode: 'copy', pattern: "${sample_id}.preprocessed_reads.fastq"

    input:
        tuple val(sample_id) , val(barcode) , val(cell_line), val(condition), val(replicate), path(fastq_file), path(locations_tsv)

    output:
        tuple val(sample_id) , val(barcode) , val(cell_line), val(condition), val(replicate), file("${sample_id}.preprocessed_reads.fastq"), path(locations_tsv)

    script:
    """

    cat ${fastq_file} | chopper --threads ${task.cpus} --quality ${params.phred_score} > ${sample_id}.preprocessed_reads.fastq

    echo -ne "\t${baseDir}/${params.outdir}/preprocessed/${sample_id}.preprocessed_reads.fastq" >> ${locations_tsv}
    """
}

process alignAndSort {
    label "aligning"

    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: "${sample_id}.bam*"

    input:
        tuple val(sample_id) , val(barcode) , val(cell_line), val(condition), val(replicate), path(fastq_file), path(locations_tsv)
        file reference_fasta
        file reference_bed

    output:
        tuple val(sample_id) , val(barcode) , val(cell_line), val(condition), val(replicate), file("${sample_id}.bam.bai"), file("${sample_id}.bam")
        path locations_tsv

    script:
    """
    minimap2 -t ${task.cpus} -ax splice  --junc-bed ${reference_bed} ${reference_fasta} ${fastq_file} | samtools sort -@ ${task.cpus} | samtools view -@ ${task.cpus} - -F 2308 -q ${params.q_score} -O BAM -o ${sample_id}.bam
    samtools index -@ ${task.cpus} ${sample_id}.bam
    echo -ne "\t${baseDir}/${params.outdir}/alignments/${sample_id}.bam\t${baseDir}/${params.outdir}/alignments/${sample_id}.bam.bai" >> ${locations_tsv}
    """
}

process quantifyIsoforms {
    label "liqa"

    publishDir "${params.outdir}/estimates", mode: 'copy', pattern: "${sample_id}.isoform_expression_estimates"

    input:
        tuple val(sample_id) , val(barcode) , val(cell_line), val(condition), val(replicate), file(bam_bai), file(bam_file)
        file locations_tsv
        path refgene_file

    output:
    	path "*.isoform_expression_estimates"
        path(locations_tsv), emit: locations_tsv

    script:
    """
    liqa -task quantify -refgene ${refgene_file} -bam ${bam_file} -out ${sample_id}.isoform_expression_estimates -max_distance 20 -f_weight 1

    echo -ne "\t${baseDir}/${params.outdir}/estimates/${sample_id}.isoform_expression_estimates" >> ${locations_tsv}
    """
}

process combiningTsvFiles {
    label "smallJob"

    publishDir "${params.outdir}/reference", mode: 'copy', pattern: "${params.phred_score}_${params.q_score}_locations.tsv"

    input:
        file locations_tsv

    output:
        file "${params.phred_score}_${params.q_score}_locations.tsv"
    
    script:
    $/
    python3 - <<EOF
    import os
    import pandas as pd

    columns = ['sample_id','cell_line','condition','replicate', 'raw_fastq', 'processed_fastq', 'bam', 'index',' isoform_expression_estimates']

    result_dataframe = pd.read_csv("${locations_tsv}",sep='\t',header=None)

    result_dataframe.columns = columns
    result_dataframe.to_csv("${params.phred_score}_${params.q_score}_locations.tsv", index=False, sep='\t')

    EOF
    /$

}

process readsAnalysis {
    executor = "local"

    publishDir "${params.outdir}/QC_metrics", mode: 'copy'

    input:
        file locations_tsv

    output:
        file "*.png"
        file "*.tsv"

    script:
    """
    ml scipy-stack
    python ${baseDir}/qc_metrics.py "${locations_tsv}" ${params.phred_score}
    """
}

process makeconditionfiles {
    executor = "local"

    publishDir "${params.outdir}/reference", mode: 'copy', pattern: "*_isoform_expression_estimates.txt"

    input:
        file locations_tsv

    output:
        path("locations.tsv"), emit: locations
        path("*_isoform_expression_estimates.txt")
    
    script:
    $/
    python3 - <<EOF
    import pandas as pd
    import sys
    
    df = pd.read_csv("${locations_tsv}",sep='\t',header=0)

    condensed_data = []
    grouped = df.groupby(['cell_line', 'condition'])

    for (cell_line, condition), group in grouped:
        estimates = group[" isoform_expression_estimates"].tolist()
        filename = f"{cell_line}_{condition}_isoform_expression_estimates.txt"
        with open(filename, 'w') as file:
            for estimate in estimates:
                file.write(estimate)
                file.write("\n")
        condensed_data.append({'cell_line': cell_line, 'condition': condition, 'filename': filename})

    condensed_df = pd.DataFrame(condensed_data)
    condensed_df.to_csv('condensed_summary.csv', index=False)
    
    unique_cell_lines = condensed_df['cell_line'].unique()
    unique_cell_lines_df = pd.DataFrame(unique_cell_lines, columns=['cell_line'])

    for index, row in unique_cell_lines_df.iterrows():
        unique_cell_lines_df.at[index,"control_filename"] = "${baseDir}" + "/" + "${params.outdir}" + f"/reference/{row.cell_line}_control_isoform_expression_estimates.txt"
        unique_cell_lines_df.at[index,"treatment_filename"] = "${baseDir}" + "/" + "${params.outdir}" + f"/reference/{row.cell_line}_treatment_isoform_expression_estimates.txt"
    
    unique_cell_lines_df.to_csv("locations.tsv",sep='\t',header=True,index=False)
    
    EOF
    /$
}

process liqaDiff {
    label "liqaDiff"

    publishDir "${params.outdir}/results", mode: 'copy', pattern: "${cell_line}_results_file" 

    input:
        tuple val(cell_line), path(control_filename), path(treatment_filename)

    output:
        tuple val(cell_line), file("${cell_line}_results_file")
    
    script:
        """
        liqa -task diff -condition_1 ${control_filename} -condition_2 ${treatment_filename} -out ${cell_line}_results_file
        """
}

Channel
    .fromPath( params.sample_sheet )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple( row.sample_id,row.barcode,row.cell_line,row.condition,row.replicate,row.fastq_dir ) }
    .set { data_ch }


workflow {

    if (params.basecalling == true) {
        fast5_ch = Channel.fromPath(params.fast5)
        fastq_ch = basecalling(fast5_ch)
    }
    else {
        fastq_ch = Channel.fromPath(params.fastq_dir)
    }

    reference_fna = file(params.reference_fna)
    reference_gtf = file(params.reference_gtf)
    reference_bed = file(params.reference_bed)

    refgene_ch = makeRefgene(reference_gtf)

    removeAdaptors(data_ch)
    preprocessingReads(removeAdaptors.out.reads)
    alignAndSort(preprocessingReads.out,reference_fna,reference_bed)
    quantifyIsoforms(alignAndSort.out,refgene_ch)
    combiningTsvFiles(quantifyIsoforms.out.locations_tsv.collectFile(name: 'sample.txt', newLine: true))
    readsAnalysis(combiningTsvFiles.out)

    // ch = Channel.fromPath("/home/bdpitica/scratch/mRNA_isoform/liqa.output/reference/10_60_locations.tsv")
    makeconditionfiles(combiningTsvFiles.out)

    makeconditionfiles.out.locations
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple( row.cell_line,row.control_filename,row.treatment_filename ) }
        .set { results_ch }
    liqaDiff(results_ch)

}
