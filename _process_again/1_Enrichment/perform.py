import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import matplotlib.pyplot as plt

# 1) load your sample‐level metadata
meta = pd.read_csv(
    '/l/users/darya.taratynova/scbb_project/_process_again/dataset_meta.csv'
).set_index('sample_id')

# 2) find all the per‐sample h5ads, pull out only microglia cells, sum to pseudobulk
input_root = '/l/users/darya.taratynova/scbb_project/_process_again'
h5ads = glob.glob(os.path.join(input_root, '**', 'ann_data_mapped', '*.h5ad'), recursive=True)

pb_list = []
for p in h5ads:
    sample = os.path.basename(p).replace('.h5ad','')
    ad = sc.read_h5ad(p)
    ad.obs['gender']    = meta.loc[sample,'gender']
    ad.obs['condition'] = meta.loc[sample,'condition']

    micro = ad[ad.obs['cell_type_marker']=='microglia']
    if micro.n_obs == 0:
        continue

    # sum counts across all microglia cells of that sample
    counts = np.array(micro.X.sum(axis=0)).ravel()
    pb_list.append(pd.Series(counts, index=micro.var_names, name=sample))

# assemble pseudobulk DataFrame
pb_df = pd.DataFrame(pb_list)
pb_df.index.name = 'sample_id'

# subset meta to just those samples
pb_meta = meta.loc[pb_df.index]

print(f"Pseudobulk matrix: {pb_df.shape[0]} samples × {pb_df.shape[1]} genes")

# 3) load your pathway definitions
gene_sets = gp.get_library(name='GO_Biological_Process_2021')
gene_set_name = 'neurogenesis (GO:0022008)'
geneset = gene_sets[gene_set_name]

import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp

def prerank_group_gsea_and_plot(counts_df, meta_df, geneset, gene_set_name,
                                group1, group2, output_prefix,
                                min_size=5, max_size=None, nperm=1000, seed=42):
    """
    Generic pairwise GSEA on pseudobulk:
      group1, group2: dicts { 'gender':…, 'condition':… }
    """
    # pull out the two sets of samples
    mask1 = (meta_df['gender']    == group1['gender'])   & \
            (meta_df['condition'] == group1['condition'])
    mask2 = (meta_df['gender']    == group2['gender'])   & \
            (meta_df['condition'] == group2['condition'])
    df1 = counts_df.loc[mask1]
    df2 = counts_df.loc[mask2]
    print(f"{group1} → {df1.shape[0]} samples,   {group2} → {df2.shape[0]} samples")

    # mean expression per gene
    m1 = df1.mean(axis=0)
    m2 = df2.mean(axis=0)

    # log2FC ranking: group1 minus group2
    logfc = np.log2(m1 + 1) - np.log2(m2 + 1)
    logfc.name = 'logFC'
    logfc.index = logfc.index.str.upper()

    # filter
    rank_series = logfc.dropna()
    avg = (m1 + m2) / 2
    keep = avg[avg >= 1].index.str.upper()
    rank_series = rank_series.loc[ rank_series.index.isin(keep) ]

    # break ties with tiny jitter
    rng   = np.random.default_rng(seed)
    noise = rng.uniform(-1e-6, 1e-6, size=rank_series.shape)
    rank_series = rank_series + noise

    # preranked GSEA
    if max_size is None:
        max_size = len(rank_series)
    pre = gp.prerank(
        rnk            = rank_series,
        gene_sets      = {gene_set_name: geneset},
        permutation_num= nperm,
        seed           = seed,
        outdir         = None,
        verbose        = False,
        min_size       = min_size,
        max_size       = max_size
    )
    res = pre.results[gene_set_name]
    print(f"→ NES = {res['nes']:.2f},  FDR = {res['fdr']:.3g}")

    # plot ES curve
    ES, hits = res['RES'], res['hits']
    fig, ax = plt.subplots(figsize=(5,3.5))
    ax.plot(ES, color='#004E89', lw=1.5)
    ax.axhline(0, color='black', ls='--', lw=0.8)
    for h in hits:
        ax.axvline(h, color='black', ymin=0.95, ymax=1.0, lw=0.6)
    ax.set(
        xlabel='Rank',
        ylabel='Enrichment score',
        title=f"{gene_set_name}\n{group1['gender']} {group1['condition']} vs\n"
              f"{group2['gender']} {group2['condition']}"
    )
    ax.title.set_fontsize(8) 
    for spine in ['top','right']:
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    fig.savefig(f"{output_prefix}.png", dpi=300)
    plt.close(fig)


# your pathway and pseudobulk from before
gene_sets      = gp.get_library(name='GO_Biological_Process_2021')
gene_set_name  = 'neurogenesis (GO:0022008)'
geneset        = gene_sets[gene_set_name]

# run the Female vs Male comparison within Disease
prerank_group_gsea_and_plot(
    counts_df     = pb_df,
    meta_df       = pb_meta,
    geneset       = geneset,
    gene_set_name = gene_set_name,
    group1        = {'gender':'Male',   'condition':'Disease'},
    group2        = {'gender':'Female', 'condition':'Disease'},
    output_prefix = 'GSEA_microglia_M_vs_F_AD'
)

# 2) Female AD vs Female Control
prerank_group_gsea_and_plot(
    counts_df     = pb_df,
    meta_df       = pb_meta,
    geneset       = geneset,
    gene_set_name = gene_set_name,
    group1        = {'gender':'Female','condition':'Disease'},
    group2        = {'gender':'Female','condition':'Control'},
    output_prefix = 'GSEA_microglia_F_AD_vs_F_Control'
)

# 3) Male AD vs Male Control
prerank_group_gsea_and_plot(
    counts_df     = pb_df,
    meta_df       = pb_meta,
    geneset       = geneset,
    gene_set_name = gene_set_name,
    group1        = {'gender':'Male','condition':'Disease'},
    group2        = {'gender':'Male','condition':'Control'},
    output_prefix = 'GSEA_microglia_M_AD_vs_M_Control'
)
 
scavenger_receptor_genes = [
    'SCARA1', 'MSR1',     'SR-AI', 'SR-AII','SCARA3', 'MSRL1',    'APC7', 'COLEC12','SCARA4',   'SRCL', 'CL-P1',
    'SCARA5', 'TESR', 'MARCO',  'SCARA2', 'SCARB1', 'SR-BI',     'CD36L1', 'CD36',   'SCARB3',    'PAS4',
    'CD68',   'SCARD1','LOX-1',  'OLR1','Dectin-1','CLEC7A','MRC1',   'CD206',     'Mannose receptor 1',
    'ASGPR1', 'CLEC4H1',   'HL-1','SCARF1', 'SREC-I', 'MEGF10', 'EMARDD','SR-PSOX','CXCL16',
    'FEEL-1', 'STAB1',     'CLEVER-1','FEEL-2', 'STAB2','CD163',  'CD163A',    'M130'
]
scavenger_gene_set = {'Scavenger_Receptors': [g.upper() for g in scavenger_receptor_genes]}

def prerank_custom_sex_gsea(counts_df, meta_df, geneset, gene_set_name,
                            condition, output_prefix):
    # 1) subset to the condition of interest
    df = counts_df.loc[meta_df['condition'] == condition]
    md = meta_df.loc[df.index]

    # 2) compute mean expression by sex
    m_f = df.loc[md['gender']=='Female'].mean(axis=0)
    m_m = df.loc[md['gender']=='Male'  ].mean(axis=0)

    # 3) log₂FC (Female − Male)
    logfc = np.log2(m_f + 1) - np.log2(m_m + 1)
    logfc.index = logfc.index.str.upper()
    rank_series = logfc.dropna()

    # 4) drop lowly expressed genes
    avg = (m_f + m_m) / 2
    keep = avg[avg >= 1].index.str.upper()
    rank_series = rank_series.loc[rank_series.index.isin(keep)]

    # 5) convert to an explicit rank (no ties!)
    #    ‘first’ gives each gene a unique integer rank
    ranked = rank_series.rank(method='first', ascending=False)
    ranked.name = 'logFC_rank'

    # 6) prerank GSEA on that unique ranking
    pre = gp.prerank(
        rnk            = ranked,
        gene_sets      = {gene_set_name: geneset},
        permutation_num=1000,
        seed           = 42,
        outdir         = None,
        verbose        = False,
        min_size       = 1,
        max_size       = len(ranked),
    )
    res = pre.results[gene_set_name]
    print(f"{condition}: Female vs Male → NES={res['nes']:.2f}, FDR={res['fdr']:.3g}")

    # 7) plot the ES curve
    ES, hits = res['RES'], res['hits']
    fig, ax = plt.subplots(figsize=(5,3.5))
    ax.plot(ES, color='#004E89', lw=1.5)
    ax.axhline(0, color='black', ls='--', lw=0.8)
    for h in hits:
        ax.axvline(h, color='black', ymin=0.95, ymax=1.0, lw=0.6)
    ax.set(
        xlabel='Rank',
        ylabel='Enrichment score',
        title=f"{gene_set_name}\nFemale vs Male in {condition}"
    )
    ax.title.set_fontsize(8) 

    for sp in ['top','right']:
        ax.spines[sp].set_visible(False)
    plt.tight_layout()
    fig.savefig(f"{output_prefix}.png", dpi=300)
    plt.close(fig)

# 3) run it for Disease
prerank_custom_sex_gsea(
    counts_df     = pb_df,
    meta_df       = pb_meta,
    geneset       = scavenger_gene_set['Scavenger_Receptors'],
    gene_set_name = 'Scavenger_Receptors',
    condition     = 'Disease',
    output_prefix = 'GSEA_scavenger_microglia_F_vs_M_Disease'
)