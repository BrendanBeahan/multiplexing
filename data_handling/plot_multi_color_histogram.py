import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os
import numpy as np

# Set the font to an available font on your system
plt.rcParams['font.family'] = 'DejaVu Sans'

# New input format including "unmapped"
# barcode sequence count
data = {}
with open('./read_counts.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 3:
            continue  # skip malformed lines
        barcode, sequence, count_str = parts
        count = int(count_str)
        
        # Store in dictionary
        if barcode not in data:
            data[barcode] = {}
        data[barcode][sequence] = data[barcode].get(sequence, 0) + count

# Extract the sorted barcodes and sequences
barcodes = sorted(data.keys(), key=lambda x: (x.lower() != 'unknown', x))  # Put unknown last if present
all_sequences = set()
for seq_dict in data.values():
    all_sequences.update(seq_dict.keys())

# Ensure unmapped reads are explicitly included
if "unmapped" not in all_sequences:
    all_sequences.add("unmapped")

all_sequences = sorted(all_sequences)

# Convert data into arrays for plotting
counts_matrix = []
for barcode in barcodes:
    seq_counts = [data[barcode].get(seq, 0) for seq in all_sequences]
    counts_matrix.append(seq_counts)

counts_matrix = np.array(counts_matrix)
total_counts = np.sum(counts_matrix)

# Choose colors for sequences, adding a distinct color for "unmapped"
colors = plt.cm.get_cmap('Paired', len(all_sequences) - 1)
sequence_colors = [colors(i) for i in range(len(all_sequences) - 1)]
sequence_colors.append("#FF5733")  # A bright orange for unmapped reads

fig, ax = plt.subplots(figsize=(8, 6))

# Bottom array for stacking
bottom = np.zeros(len(barcodes))
bars = []

# Plot each sequence as a segment in the stacked bar
for i, seq in enumerate(all_sequences):
    seq_counts = counts_matrix[:, i]
    bar = ax.bar(barcodes, seq_counts, bottom=bottom, color=sequence_colors[i], width=0.6, label=seq)
    bars.append(bar)
    bottom += seq_counts

# Y-axis formatting
ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{int(x):,}'))

# Label axes
ax.set_xlabel('Barcode', fontsize=14)
ax.set_ylabel('Number of Reads', fontsize=14)

# Customize ticks
plt.xticks(rotation=0, fontsize=12)
plt.yticks(fontsize=12)

# Removing gridlines and unnecessary spines
ax.grid(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add data labels on top of the stacked bar for total counts
for i, barcode in enumerate(barcodes):
    yval = sum(counts_matrix[i])
    percentage = yval / total_counts * 100 if total_counts > 0 else 0
    ax.text(
        i, 
        yval + max(np.sum(counts_matrix, axis=1)) * 0.01,
        f'{int(yval):,}\n({percentage:.1f}%)',
        ha='center',
        va='bottom',
        fontsize=12
    )

# Adjust y-axis limits to provide space for labels
ax.set_ylim(0, max(np.sum(counts_matrix, axis=1)) * 1.15)

# Add a legend outside the plot area to the right
ax.legend(
    title='Sequence',
    fontsize=12,
    title_fontsize=12,
    loc='upper center',  # Centered at the bottom
    bbox_to_anchor=(0.5, -0.15),  # Below the plot
    ncol=len(all_sequences)  # Display legend items in a single row
)
plt.tight_layout(rect=[0, 0.15, 1, 1])  # Adjust layout to make space below the plot

plt.savefig('./read_counts_histogram.png', dpi=300, bbox_inches='tight')
plt.show()
