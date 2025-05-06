import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1) point these at your four CSVs
files = {
    'F_AD_vs_F_Control': '/l/users/darya.taratynova/scbb_project/_process_again/1_DE/F_AD_vs_F_Control.csv',
    'M_AD_vs_M_Control': '/l/users/darya.taratynova/scbb_project/_process_again/1_DE/M_AD_vs_M_Control.csv',
    'F_vs_M_Control':    '/l/users/darya.taratynova/scbb_project/_process_again/1_DE/F_vs_M_Control.csv',
    'F_vs_M_Disease':    '/l/users/darya.taratynova/scbb_project/_process_again/1_DE/F_vs_M_Disease.csv',
}

# 2) load & pick significant genes (padj<0.05)
sig = {}
for name, path in files.items():
    df = pd.read_csv(path)
    up   = set(df.loc[(df.padj < 0.05) & (df.log2FoldChange > 0), 'gene'])
    down = set(df.loc[(df.padj < 0.05) & (df.log2FoldChange < 0), 'gene'])
    sig[name] = {'Up': up, 'Down': down}

# 3) build a "long" DataFrame for seaborn
rows = []
for cmp, sets in sig.items():
    rows.append({'cmp': cmp, 'Direction': 'Up',   'Count': len(sets['Up'])})
    rows.append({'cmp': cmp, 'Direction': 'Down', 'Count': len(sets['Down'])})
bar_df = pd.DataFrame(rows)

# 4) rename to short labels & force the order
label_map = {
    'F_AD_vs_F_Control': 'AD-F',
    'M_AD_vs_M_Control': 'AD-M',
    'F_vs_M_Control':    'Ctrl-F',
    'F_vs_M_Disease':    'Ctrl-M',
}
order = ['AD-F','Ctrl-F','Ctrl-M','AD-M']
bar_df['Label'] = bar_df['cmp'].map(label_map)
bar_df['Label'] = pd.Categorical(bar_df['Label'],
                                 categories=order,
                                 ordered=True)

# 5) plot
fig, ax = plt.subplots(figsize=(6,4))
sns.barplot(
    data=bar_df,
    x='Label',
    y='Count',
    hue='Direction',
    palette={'Down':'#4C72B0','Up':'#DD8452'},
    dodge=True,
    ax=ax
)

# 6) annotate counts on top
ymax = bar_df.Count.max()
for p in ax.patches:
    h = p.get_height()
    if h>0:
        ax.text(
            p.get_x() + p.get_width()/2,
            h + ymax*0.01,
            f"{int(h)}",
            ha='center',
            va='bottom',
            fontsize= 9,  # adjust size if you like
        )

# 7) cosmetics
ax.set_ylabel("Number of DE genes")
ax.set_xlabel("")            # no x-axis title
ax.legend(title="",         # no legend title
          loc="upper right",
          frameon=False)
plt.tight_layout()
plt.show()
plt.close(fig)
fig.savefig("DE_genes_barplot.png", dpi=300)