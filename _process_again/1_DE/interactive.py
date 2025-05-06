import pandas as pd
import numpy as np
import plotly.graph_objects as go

# 1) Load your DE results
path = '/l/users/darya.taratynova/scbb_project/_process_again/1_DE/results/DE_gender_within_control.csv'
de = pd.read_csv(path, index_col=0)
de['-log10_padj'] = -np.log10(de['padj'].replace(0, np.nan))

# 2) Define thresholds and split
fc_thresh, padj_thresh = 1.0, 0.05
up   = de[(de['padj']<padj_thresh) & (de['log2FoldChange']>= fc_thresh)]
down = de[(de['padj']<padj_thresh) & (de['log2FoldChange']<=-fc_thresh)]
ns   = de.drop(up.index).drop(down.index)

# 3) Build the figure
fig = go.Figure()

# Not significant
fig.add_trace(go.Scatter(
    x=ns['log2FoldChange'], y=ns['-log10_padj'],
    mode='markers',
    name='NotSig',
    marker=dict(color='lightgray', size=6),
    hovertext=ns.index,
    hoverinfo='text'
))

# Up-regulated (with text labels)
fig.add_trace(go.Scatter(
    x=up['log2FoldChange'], y=up['-log10_padj'],
    mode='markers+text',
    text=up.index,
    textposition='top center',
    textfont=dict(size=8, color='red'),
    name='Up (padj<0.05 & log2FC≥1)',
    marker=dict(color='red', size=6),
    hovertext=up.index,
    hoverinfo='text'
))

# Down-regulated (with text labels)
fig.add_trace(go.Scatter(
    x=down['log2FoldChange'], y=down['-log10_padj'],
    mode='markers+text',
    text=down.index,
    textposition='bottom center',
    textfont=dict(size=8, color='blue'),
    name='Down (padj<0.05 & log2FC≤-1)',
    marker=dict(color='blue', size=6),
    hovertext=down.index,
    hoverinfo='text'
))

# 4) Add threshold lines
fig.add_hline(y=-np.log10(padj_thresh), line_dash='dash', line_color='black')
fig.add_vline(x= fc_thresh,                line_dash='dash', line_color='black')
fig.add_vline(x=-fc_thresh,                line_dash='dash', line_color='black')

# 5) Layout tweaks
fig.update_layout(
    title='Volcano Plot of DESeq2 Results',
    xaxis_title='Log₂ Fold Change',
    yaxis_title='-log₁₀(padj)',
    legend=dict(itemsizing='constant')
)

# 6) Save & show
html_path = '/l/users/darya.taratynova/scbb_project/_process_again/1_DE/results/volcano_interactive_gender_control.html'
fig.write_html(html_path)
print("→ Saved interactive volcano to:", html_path)

fig.show()
