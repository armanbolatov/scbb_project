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

    micro = ad[ad.obs['cell_type_marker']=='astocytes']
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
pathway_names = [
    'cellular response to chemical stress (GO:0062197)',
    'response to zinc ion (GO:0010043)',
    'positive regulation of interleukin-8 production (GO:0032757)',
    'negative regulation of inclusion body assembly (GO:0090084)',       # e.g. “signalling” vs “signaling”
    'regulation of nucleotide-binding oligomerization domain containing 2 signaling pathway (GO:0070432)',           # maybe plural “processes”
    'negative regulation of transcription, DNA-templated (GO:0045892)',
    'regulation of angiogenesis (GO:0045765)',
    'regulation of trans-synaptic signaling (GO:0099177)',
    'fructose 6-phosphate metabolic process (GO:0006002)',
    'regulation of defense response (GO:0031347)',
    'regulation of response to biotic stimulus (GO:0002831)'
]
# for term in gene_sets:
#     if 'chemical' in term.lower():
#         print(term)

#     if 'transcription' in term.lower():
#         print(term)
# for go_id in ['GO:0099536','GO:0099537','GO:0098609',
#               'GO:0016043','GO:0008152','GO:0046034']:
#     matches = [term for term in gene_sets if go_id in term]
#     print(go_id, '→', matches or '⚠️ not present at all')

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
    ax.plot(ES, color='#FF6B35', lw=1.5)
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


# 3) Option A: pairwise Male vs Female in Disease
for pname in pathway_names:
    if pname not in gene_sets:
        print(f"⚠️ pathway not found in library: {pname}")
        continue

    geneset = gene_sets[pname]
    prerank_group_gsea_and_plot(
        counts_df     = pb_df,
        meta_df       = pb_meta,
        geneset       = geneset,
        gene_set_name = pname,
        group1        = {'gender':'Male',   'condition':'Disease'},
        group2        = {'gender':'Female', 'condition':'Disease'},
        output_prefix = f"_astocytes/{pname}_GSEA_M_vs_F_{pname.replace(' ','_')}"
    )

# 4) Option B: Female AD vs Female Control (and likewise for males)
for sex in ['Female','Male']:
    for pname in pathway_names:
        if pname not in gene_sets:
            continue
        prerank_group_gsea_and_plot(
            counts_df     = pb_df,
            meta_df       = pb_meta,
            geneset       = gene_sets[pname],
            gene_set_name = pname,
            group1        = {'gender':sex, 'condition':'Disease'},
            group2        = {'gender':sex, 'condition':'Control'},
            output_prefix = f"_astocytes/{pname}_GSEA_{sex}_AD_vs_{sex}_Control_{pname.replace(' ','_')}"
        )