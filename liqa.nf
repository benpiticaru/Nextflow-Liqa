#!/usr/bin/env nextflow

params.fast5 = "/home/bdpitica/scratch/mRNA_isoform/minion_combined/pod5"
params.sample_sheet = "/home/bdpitica/scratch/mRNA_isoform/reference/liqa_reference.tsv"
params.phred_score = 10
params.q_score = 60
params.outdir = "final_liqa.output"
params.basecalling = false
params.fastq_dir = "/home/bdpitica/scratch/mRNA_isoform/guppy.combined.out/pass"
params.reference_gtf = "/home/bdpitica/scratch/mRNA_isoform/reference/Homo_sapiens_GRCh38_110.gtf"
params.reference_fna = "/home/bdpitica/scratch/mRNA_isoform/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.reference_bed = "/home/bdpitica/scratch/mRNA_isoform/reference/Homo_sapiens_GRCh38_110.bed"

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

process goAnalysis {
    label "smallJob"

    publishDir "${params.outdir}/Go_analysis", mode: 'copy' 

    input:
        tuple val(cell_line), file(results_file)

    output:

    
    script:
        $/
        python3 - <<EOF
        import pandas as pd
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import numpy as np
        import seaborn as sns
        import textwrap
        from statsmodels.stats import multitest
        from gprofiler import GProfiler
        import matplotlib.ticker as ticker

        def analyze_and_plot_gene_data(results_file_path, cell_line, significance_level=0.05):
            df = pd.read_csv(results_file_path, sep="\t", header=None, names=["gene_name", "p-value"])
            _, q_values, _, _ = multitest.multipletests(df['p-value'], alpha=significance_level, method='fdr_bh')
            df['q-value'] = q_values
            df = df.sort_values(by="p-value").query("`q-value` < @significance_level").sort_values(by="q-value")
            df.to_csv(f"{cell_line}_adj_p-values.tsv", sep='\t', index=False)
            
            gp = GProfiler(return_dataframe=True)
            go_terms_df = gp.profile(organism='hsapiens', query=df["gene_name"].tolist(), sources=["GO:MF", "GO:CC", "GO:BP", "KEGG"])
            go_terms_df.to_csv(f"{cell_line}_go_terms.tsv", sep='\t')
            
            cmap = mpl.cm.bwr_r
            norm = mpl.colors.Normalize(vmin=0, vmax=significance_level)
            mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
            palette = mapper.to_rgba(go_terms_df['p_value']).tolist()

            fig, ax = plt.subplots(figsize=(6, 4))
            barplot = sns.barplot(data=go_terms_df, x='recall', y='name', hue='name', palette=palette, ax=ax, edgecolor='black', linewidth=1, legend=False)

            ax.set_yticks(np.arange(len(go_terms_df['name'])))
            ax.set_yticklabels([textwrap.fill(e, 26) for e in go_terms_df['name']], fontsize=12)
            ax.yaxis.set_major_locator(ticker.FixedLocator(ax.get_yticks()))
            ax.set_xlabel('Recall', fontsize=12)
            ax.set_ylabel('Name', fontsize=12)

            cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])
            cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')
            cbar.set_ticks(np.linspace(0, significance_level, num=5))  
            cbar.set_ticklabels(['{:.3f}'.format(v) for v in np.linspace(0, significance_level, num=5)])
            cbar.set_label('adj. p-value', fontsize=12)

            plt.savefig(f"{cell_line}_GO_analysis.png", format="png", bbox_inches="tight")
            plt.close()

        analyze_and_plot_gene_data(${results_file},${cell_line})

        EOF
        /$
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

    makeconditionfiles(combiningTsvFiles.out)
    makeconditionfiles.out.locations
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple( row.cell_line,row.control_filename,row.treatment_filename ) }
        .set { results_ch }
    liqaDiff(results_ch)
    goAnalysis(liqaDiff.out)
}