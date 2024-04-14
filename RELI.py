"""
A basic reimplementation of the RELI (Regulatory Element Locus Interection)
algorithm as created by the Weirach Lab: https://www.nature.com/articles/s41588-018-0102-3
"""

"""
Required inputs:
- ChIP-seq index file - An overview of the info for each ChIP-seq sample analyzed (ex. name, lab, cell, label, etc).
- LD blocks - This contains a list of 
- Target - The ChIP-seq sample to observe.
- Build - The list of where each chromosome starts for the genome being analyzed.
- Null file - 
- DBSNP - A database containing a list of all commmon SNPs for the genome.
- rep - Number of times to repeat the algorithm
- 
"""
class RELI:
    def read_ld(self, src):
        # Read this file line-by-line
        f = open(src, "r+")
        for line in f.readlines():
            # For each line, take the first value
            pass
    def __init__(
        self,
        num_reps,
        ld_file
    ):
        # Read and initialize all the information given
        self.num_reps = num_reps
        # Read the LD file given
        self.read_ld()
    def load_ld(self):
        """load LD SNPs, pushing them into ld_vectors

        :param p1: describe about parameter p1
        :param p2: describe about parameter p2
        :param p3: describe about parameter p3
        :return: describe what it returns
        """
        
    def simulate(self):
        """does blah blah blah.

        :param p1: describe about parameter p1
        :param p2: describe about parameter p2
        :param p3: describe about parameter p3
        :return: describe what it returns
        """
        # (1) The first time around, iterate over all LD vectors,
        for ld_vec in self.ld_vectors:
            pass
        # (2) Iterate repmax times over the sequence
        for i in range(self.num_reps):
            pass