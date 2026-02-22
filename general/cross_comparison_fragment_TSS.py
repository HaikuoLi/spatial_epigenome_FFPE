import glob
import gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pysam

# Find all fragments.bed.gz files in the current directory
bed_files = glob.glob("*.bed.gz")

# Define custom colors for the samples (modify as needed)
custom_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf','black']


# Ensure the number of colors matches the number of samples
if len(custom_colors) < len(bed_files):
    raise ValueError("Please provide enough colors for all samples.")

# Initialize a list to store fragment sizes for each sample
all_data = []

# Process each BED file
for bed_file in bed_files:
    sample_name = bed_file.replace(".bed.gz", "")  # Extract sample name
    fragment_sizes = []
    
    # Read the gzipped BED file
    with gzip.open(bed_file, "rt") as f:
        for line in f:
            if line.startswith("#"):  # Skip header lines
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 3:  # Ensure the line has at least chr, start, end
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                    fragment_size = end - start
                    if fragment_size > 0:  # Only include positive fragment sizes
                        fragment_sizes.append(fragment_size)
                except ValueError:
                    continue  # Skip lines with invalid numeric values
    
    # Create a DataFrame for the current sample
    df = pd.DataFrame({"Fragment Size": fragment_sizes, "Sample": sample_name})
    all_data.append(df)

# Combine all sample data into a single DataFrame
combined_df = pd.concat(all_data, ignore_index=True)

order = [
    'Fresh mouse E11', 'Fresh human tonsil', 
    'FFPE mouse brain', 'FFPE mouse colon', 'FFPE human LN', 
    'Fresh mouse E11 (H3K27me3)', 'Fresh mouse E11 (H3K4me3)',
    'FFPE human LN (H3K27me3)-1', 'FFPE human LN (H3K27me3)-2', 
    'FFPE human LN (H3K4me3)'
]
combined_df['Sample'] = pd.Categorical(combined_df['Sample'], categories=order, ordered=True)
combined_df = combined_df.sort_values('Sample')

plt.figure(figsize=(4,4))
ax = sns.histplot(
    data=combined_df,
    x="Fragment Size",
    hue="Sample",
    palette=custom_colors[:len(bed_files)],
    element="poly",
    stat="count",
    common_norm=False,
    bins=400,
    log_scale=(False, False),
    fill=False,  # No fill, only outline
    alpha=0.85
)

ax.legend_.remove()  # same color scheme as in TSSplot

plt.xlabel("Fragment Insert Size (bp)", fontsize=10.5)
plt.ylabel("Count", fontsize=11)
plt.xlim(0, 500)
plt.tight_layout()
plt.savefig('fragment_size_distribution.png', dpi=300)



##########################



import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def plot_tss_enrichment_multiple(distances_list, labels, window=2000, bin_size=10, background_start=1500, background_end=2000, sigma=3, ymax=1.0, output_file="tss_enrichment.pdf"):
    """Plot TSS enrichment profiles with count normalization, background normalization, smoothing, and min-max normalization."""
    if len(distances_list) != len(labels):
        raise ValueError("Number of distance arrays must match number of labels")
    
    # Create bins for histogram
    bins = np.arange(-window, window + bin_size, bin_size)
    plt.figure(figsize=(10, 6))
    
    for distances, label in zip(distances_list, labels):
        if len(distances) == 0:
            print(f"Warning: Empty distance array for {label}")
            continue
        
        # Calculate histogram, normalized by total counts
        hist, bin_edges = np.histogram(distances, bins=bins, density=True)
        hist = hist * len(distances) / (bin_size * len(distances))  # Adjust for bin size and total counts
        
        # Background normalization
        background_bins = (bin_edges[:-1] >= background_start) | (bin_edges[:-1] <= -background_start)
        background_counts = hist[background_bins]
        background_mean = np.mean(background_counts) if len(background_counts) > 0 else 1.0
        hist_normalized = hist / background_mean if background_mean > 0 else hist
        
        # Apply Gaussian smoothing
        hist_smoothed = gaussian_filter1d(hist_normalized, sigma=sigma)

        # Normalize to ymax and min=0
        max_height = np.max(hist_smoothed) if len(hist_smoothed) > 0 else 1.0
        min_height = np.min(hist_smoothed) if len(hist_smoothed) > 0 else 0.0
        if max_height > min_height:
            # Scale to range [0, ymax]
            hist_smoothed = (hist_smoothed - min_height) * (ymax / (max_height - min_height))
        else:
            print(f"Warning: Max height ({max_height:.4f}) <= min height ({min_height:.4f}) for {label}, using unscaled data")

        plt.plot(bin_edges[:-1], hist_smoothed, label=label, linewidth=1.8)  # Increased line thickness
    
    # Customize plot
    plt.xlabel('Distance to TSS (bp)', fontsize=13)
    plt.ylabel('Normalized Fragment Inserts', fontsize=13.5)
    plt.title('')
    plt.grid(True)
    plt.axvline(0, color='red', linestyle='--', alpha=0.5)
    #plt.legend(bbox_to_anchor=(1.05, 0.75), loc='upper left', borderaxespad=0., fontsize=11.5)  # Legend outside right
    legend = plt.legend(bbox_to_anchor=(1.05, 0.8), loc='upper left', borderaxespad=0., fontsize=11.5, frameon=False)
    for line in legend.get_lines():
        line.set_linewidth(4)  # Set legend line thickness
    plt.tight_layout()
    
    # Save as high-resolution PDF
    plt.savefig(output_file, format='png', dpi=300, bbox_inches='tight')
    plt.show()



def load_fragments(fragment_file):
    """Load fragment.bed or fragment.bed.gz file, skip lines starting with #, and extract fragment ends."""
    # Check if file is gzipped based on extension
    is_gzipped = fragment_file.endswith('.gz')
    
    # Read BED file (tab-separated, no header, skip lines starting with #)
    if is_gzipped:
        with gzip.open(fragment_file, 'rt') as f:
            fragments = pd.read_csv(f, sep='\t', header=None, 
                                   names=['chr', 'start', 'end', 'cell_id', 'count'],
                                   comment='#')
    else:
        fragments = pd.read_csv(fragment_file, sep='\t', header=None, 
                               names=['chr', 'start', 'end', 'cell_id', 'count'],
                               comment='#')
    
    # Extract fragment ends (start and end positions)
    ends = []
    for _, row in fragments.iterrows():
        ends.append({'chr': row['chr'], 'pos': row['start']})
        ends.append({'chr': row['chr'], 'pos': row['end']})
    return pd.DataFrame(ends)

from scipy.ndimage import gaussian_filter1d
def calculate_distances(fragments, tss, window=2000):
    """Calculate distance of fragment ends to nearest TSS."""
    distances = []
    for chr_name in fragments['chr'].unique():
        frag_chr = fragments[fragments['chr'] == chr_name]['pos'].values
        tss_chr = tss[tss['chr'] == chr_name]['tss_pos'].values
        
        if len(tss_chr) == 0 or len(frag_chr) == 0:
            continue
            
        for pos in frag_chr:
            dists = pos - tss_chr
            dists = dists[np.abs(dists) <= window]
            if len(dists) > 0:
                min_dist = dists[np.argmin(np.abs(dists))]
                distances.append(min_dist)

def shareseq_load_tss(tss_file):
    tss = pd.read_csv(tss_file, sep='\t', header=None, 
                      names=['chr', 'tss_pos','tss_pos2' ,'strand'])
    return tss
    
#one example of how to calculate TSS distance from a bed file
shareseq_tss_mm10 = shareseq_load_tss('./SHAREseq_tss/mm10.TSS.bed')
fragments_fresh_ME11 = load_fragments('fresh_GSM5238385_ME11_50um.downsampled.fragment.bed.gz')
distances_fresh_ME11 = calculate_distances(fragments_fresh_ME11, shareseq_tss_mm10, window=2000)


# Example usage
plot_tss_enrichment_multiple(
    distances_list=[distances_fresh_ME11, distances_fresh_htonsil,
                    distances, distances_mouse_colon, distances_malt, 
                    distances_fresh_ME11_H3K27me3,distances_fresh_ME11_H3K4me3,
                    distances_H3K27me3_1, distances_H3K27me3_2, distances_H3K4me3_1],
    labels=['Fresh mouse E11', 'Fresh human tonsil', 
            'FFPE mouse brain', 'FFPE mouse colon', 'FFPE human LN', 
            'Fresh mouse E11 (H3K27me3)','Fresh mouse E11 (H3K4me3)',
            'FFPE human LN (H3K27me3)-1', 'FFPE human LN (H3K27me3)-2', 'FFPE human LN (H3K4me3)'],
    window=2000,
    bin_size=10,
    background_start=1500,
    background_end=2000,
    sigma=6,
    ymax=1.0,
    output_file="20250705_tss_enrichment_plot.png"
)