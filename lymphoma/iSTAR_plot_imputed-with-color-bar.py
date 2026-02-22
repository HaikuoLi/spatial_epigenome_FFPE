import sys

import numpy as np
import matplotlib.pyplot as plt

from utils import load_pickle, save_image, read_lines, load_image
# from visual import cmap_turbo_truncated


def plot_super(
        x, outfile, gene_name, underground=None, truncate=None):
    ##    x, outfile, underground=None, truncate=None):
    
    x = x.copy()
    mask = np.isfinite(x)

    if truncate is not None:
        x -= np.nanmean(x)
        x /= np.nanstd(x) + 1e-12
        x = np.clip(x, truncate[0], truncate[1])

    x -= np.nanmin(x)
    x /= np.nanmax(x) + 1e-12

    cmap = plt.get_cmap('turbo')
    # cmap = cmap_turbo_truncated
    if underground is not None:
        under = underground.mean(-1, keepdims=True)
        under -= under.min()
        under /= under.max() + 1e-12

    img = cmap(x)[..., :3]
    if underground is not None:
        img = img * 0.5 + under * 0.5
    img[~mask] = 1.0
    img = (img * 255).astype(np.uint8)
    ## save_image(img, outfile)

    ## the below are added/modified by Mingyu Yang(09/11/2024)

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(12, 12))  # Adjust figure size as needed
    im = ax.imshow(img)

    # Add title with gene name
    ax.set_title(gene_name, fontsize=24, pad=20)

    # Hide axis ticks and borders
    ax.set_xticks([])
    ax.set_yticks([])
    plt.axis('off')  # Turn off axis and border

    # Add colorbar to the right of the image with shrink factor to control size
    cbar = plt.colorbar(im, ax=ax, orientation='vertical', fraction=0.046, pad=0.04, shrink=0.6)

    # Save the image without any border
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight', pad_inches=0)
    plt.close(fig)

def main():

    prefix = sys.argv[1]  # e.g. 'data/her2st/B1/'
    gene_names = read_lines(f'{prefix}gene-names.txt')
    mask = load_image(f'{prefix}mask-small.png') > 0

    for gn in gene_names:
        cnts = load_pickle(f'{prefix}cnts-super/{gn}.pickle')
        cnts[~mask] = np.nan
        # Generate the output PNG file path
        output_file = f'{prefix}cnts-super-plots/{gn}.png'
        plot_super(cnts, output_file, gene_name=gn)
        ## plot_super(cnts, f'{prefix}cnts-super-plots/{gn}.png')

if __name__ == '__main__':
    main()
