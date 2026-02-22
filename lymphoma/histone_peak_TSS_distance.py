def read_bed(file_path, is_tss=False):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if is_tss and len(fields) >= 4:
                data.append([fields[0], int(fields[1]), int(fields[2]), fields[3]])  # chrom, start, end, strand
            elif not is_tss and len(fields) >= 10:  # Assuming narrowPeak format
                data.append([fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4]), fields[5], float(fields[6])])
    return pd.DataFrame(data, columns=['chrom', 'start', 'end', 'strand'] if is_tss else ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue'])


# Function to calculate distance from peak center to nearest TSS considering strand
def find_nearest_tss(peak_start, peak_end, peak_chrom, tss_df):
    peak_center = (peak_start + peak_end) / 2
    same_chrom = tss_df[tss_df['chrom'] == peak_chrom]
    if same_chrom.empty:
        return None
    min_dist = float('inf')
    for idx, tss in same_chrom.iterrows():
        tss_start = tss['start']
        strand = tss['strand']
        if strand == '+':
            dist = tss_start - peak_center  # Distance from peak center to TSS
        else:  # strand == '-'
            dist = peak_center - tss_start  # Distance from peak center to TSS (reversed for negative strand)
        if abs(dist) < abs(min_dist):
            min_dist = dist
    return min_dist

narrowpeak_df2 = read_bed("../CT3_0602_H3K4me3/signac_output/MACS2/CT3_0602_peaks.narrowPeak", is_tss=False)
# Calculate distances and signal values for all peaks
distances2 = []
signal_values2 = []
for idx, peak in narrowpeak_df2.iterrows():
    dist = find_nearest_tss(peak['start'], peak['end'], peak['chrom'], tss_df)
    if dist is not None:
        distances2.append(dist)
        signal_values2.append(peak['signalValue'])  # Using signalValue as intensity metric

plt.figure(figsize=(4,6))
plt.scatter(distances2[:6000], range(len(distances2[:6000])), c=signal_values2[:6000], cmap='Purples', s=3, vmin=0, vmax=2)
cbar = plt.colorbar(label='Signal Value', fraction=0.03, pad=0.02)  # Shrink colorbar
cbar.set_ticks([0, 0.5, 1, 1.5, 2])  # Set tick positions
cbar.set_ticklabels(['0', '0.25', '0.5', '0.75', '1'])  # Custom labels scaled to 0-1
plt.grid(False)  # Remove grid
plt.xlabel('Distance from TSS (bp)')
plt.ylabel('')
plt.title('H3K4me3 peaks (DLBCL)')
plt.ylim(len(distances2[:6000]), 0)
plt.xlim(-10000, 10000)
plt.tight_layout()
plt.savefig('narrowpeakTSS.png', dpi=300)
plt.show()

np.save('TSSdistance/distances_H3K4me3.npy', distances2)
np.save('TSSdistance/signal_values_H3K4me3.npy', signal_values2)