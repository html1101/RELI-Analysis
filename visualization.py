# %% [markdown]
# # Graphs for Paper

# %% [markdown]
# ## Pt. 1: Visualizing ChIP-seq Data
# We have a lot of ChIP-seq data. Let's take a look at what exactly it looks like.

# %%
import numpy as np
import RELI as R
from RELI import RELI, LoadedData
R.DEBUG = False

# %%
# Get list of all the locations of these ChIP-seqs + plot all this
# sample = "hg19_0697"

# %%
# plt.hist(new_data, bins=num_bins)
# plt.savefig(F"{sample}.png")
# plt.show()

# %%

# Now each of these bins will be handed off to a thread, which will deal with loading the data
def handle_bin(data, bin, thread_id):
    for target in bin:
        print(F"Thread {thread_id} running on target {target}")
        # Run RELI - we're saving the results into files so no saving happens here
        RELI(data, target)

if __name__ == '__main__':
    from multiprocessing import Process
    print("Main line starting multithreading processes")
    snp_file = "mas/type3/Vitiligo.snp"
    # We first load in the data we use for all the ChIP-seq files before beginning analysis.
    data = LoadedData("SLE",
        snp_file,
        1000,
        # No LD file supplying
        None,
        "sample_data/ChIPseq.index",
        given_species = "sample_data/GenomeBuild/hg19.txt",
        output_dir="mas_type3")

    # Iterate over all the ChIP-seq options loaded in
    # Number of threads we're going to use for processing
    num_bins = 10

    # Place these into 5 bins, where we will be performing multithreading
    bins = []
    chip_values = list(data.chip_seq_index.keys())
    count = int(len(chip_values) / (num_bins - 1)) if num_bins > 1 else len(chip_values)
    while len(chip_values) > 0:
        bins.append(chip_values[:count])
        del chip_values[:count]

    threads = []
    necessary_info = data.necessary_info()
    for i, bin in enumerate(bins):
        t = Process(target=handle_bin, args=[necessary_info, bin, i])
        threads.append(t)
        t.start()
    
    # Now join all the threads to the main thread
    for thread in threads:
        thread.join()
    
    print("We're done, yay!")

