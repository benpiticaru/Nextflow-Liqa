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


# Example usage
results_file_path = "/home/bdpitica/scratch/mRNA_isoform/liqa.output/results/MCF7_results_file"
cell_line = "MCF7"
analyze_and_plot_gene_data(results_file_path, cell_line)
