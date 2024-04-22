"""
A basic reimplementation of the RELI (Regulatory Element Locus Interection)
algorithm as created by the Weirach Lab: https://www.nature.com/articles/s41588-018-0102-3
"""
# Imports
# Colored printing
from colorama import Fore, Back, Style
import sys, os, time, math
from time import sleep
import numpy as np
import random, scipy
from operator import itemgetter, attrgetter
from scipy import special as sp

def LINE():
    return sys._getframe(1).f_lineno

def quitting(line: int):
    """Quit the program with failure
    :param line: Line error occurred on
    """
    print(F"Error occurred on {line}!\nQuitting...")
    # Exit with nonzero status code
    exit(-1)

# Type of statistical model used
[NORMAL, EMPIRICAL, PHASETYPE, BINOMIAL, DEFAULT] = range(0, 5)

# Sort SNPs consecutively, first by chromosome, then by snp_start
def snpSort(lhs, rhs):
	return (lhs["snp_chr"] < rhs["snp_chr"]) or ((lhs["snp_chr"] == rhs["snp_chr"]) and (lhs["snp_start"] < rhs["snp_start"]))

# Source: https://stackoverflow.com/questions/56075838/how-to-generate-the-values-for-the-q-function
def qfunc(x):
    return 0.5-0.5*sp.erf(x/math.sqrt(2))

def overlap_exists(bed, snp):
    return ((bed["bed_start"] <= snp["snp_start"] and bed["bed_end"] >= snp["snp_end"]) or
    (bed["bed_start"] >= snp["snp_start"] and bed["bed_start"] + 1 < snp["snp_end"]) or
    (bed["bed_end"] > snp["snp_start"] + 1 and bed["bed_end"] <= snp["snp_end"]) or
    (bed["bed_start"] >= snp["snp_start"] and bed["bed_end"] <= snp["snp_end"]))

# ChIP-seq file information
class BEDSample:
    # Read file given its PATH
    def readFile(self, src: str):
        sys.stdout.write(F"\r[   ] Reading BED file (file: {src})     ")
        sys.stdout.flush()
        # (1) Read file
        for line in open(src):
            # TODO: use binding strength?
            [bed_chr, bed_start, bed_end] = line.replace("\n", "").replace("\n", "").split("\t")[:3]
            bed_start = int(bed_start)
            bed_end = int(bed_end)

            # Take the length by end - start
            length = bed_end - bed_start
            entry = {
                "bed_start": bed_start,
                "bed_end": bed_end,
                "bed_chr": bed_chr,
                "length": length
            }
            # Push the list of lengths
            self.lengths.append(length)
            # Push this into a list of all the data we know about
            self.data.append(entry)
        
        """
        if (inVal){
            vector<RELI::LD> tLDVec;
            for (auto k : this->myData){
                RELI::LD tld;
                RELI::SNP t;
                t.snp_chr = k.bed_chr;
                t.snp_start = k.bed_start;
                t.snp_end = k.bed_end;
                t.snp_name = "na";
                tld.keySNP = t;
                tld.mySNP.push_back(t);
                tld.dis2keySNP.push_back(0);

                tLDVec.push_back(tld);
            }
            std::default_random_engine tSeed(std::chrono::system_clock::now().time_since_epoch().count()); //RNG seed 
            std::uniform_int_distribution<unsigned int> tGen(0, (RELI::bg_null_model_data.bin0.size() - 1));// RNG generator
            for (auto k = tLDVec.begin(); k != tLDVec.end(); ++k){
                bool tGood;
                unsigned int tIndex;
                RELI::SNP tKeySNP;
                tKeySNP.length = k->keySNP.length;

                while (tGood != true){
                    tIndex = tGen(tSeed);
                    tGood = RELI::SNPfit(*k, tKeySNP, RELI::bg_null_model_data.bin0.at(tIndex),
                        RELI::chromosome_strucuture, RELI::chromosome_strucuture_val);	// fit the new key SNP detail info into simulated location
                }
                bed3col t;
                t.bed_chr = tKeySNP.snp_chr;
                t.bed_start = tKeySNP.snp_start;
                t.bed_end = tKeySNP.snp_end;

                this->myData_bgnull.push_back(t);
                this->myData = this->myData_bgnull;  
                this->myData_bgnull.clear();
            }
        }
        """
        # sort(this->myData.begin(), this->myData.end());
        # Sort according to length
        # arg_ordering = np.argsort(self.lengths)
        self.data.sort(key=lambda x: (x["bed_chr"], x["bed_start"]))
        self.lengths.sort()
        # self.lengths = np.array(self.lengths)[arg_ordering]
        self.median_length = self.lengths[int(len(self.lengths) / 2)]
        
        print(F"\r[ {Fore.GREEN + 'x' + Style.RESET_ALL} ] Parsed BED file (file: {src})")
    def makeIndex(self):
        # At the loc of the first target file entry set ind to 0
        bed_chr_name = self.data[0]["bed_chr"]
        self.index[bed_chr_name] = 0
        prev_val = bed_chr_name
        # Iterate through + find distance from first entry
        for (ind, entry) in enumerate(self.data):
            if entry["bed_chr"] != prev_val:
                # Distance rn is simply abs value
                prev_val = entry["bed_chr"]
                self.index[prev_val] = ind
    # Initialize given a ChIP-seq input file
    def __init__(self, src: str):
        # List of 
        self.index = {}
        self.lengths = []
        self.data = []
        # Read file
        self.readFile(src)
        # Make index
        self.makeIndex()

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
        sys.stdout.write(F"\r[   ] Creating species map (file: {src})")
        sys.stdout.flush()
        if not src:
            print(F"\r[ {Fore.RED + 'x' + Style.RESET_ALL} ] Failure reading LD file; no LD file supplied")
            quitting(LINE())
        # Read this file line-by-line
        for line in open(src, "r+"):
            # For each line, convert to LD template, which
            # contains listSNP and keySNP
            # listSNP contains a list of all SNPs associated with this particular SNP.
            if len(line) > 0:
                key_snp, corresponding_snps = line.replace("\n", "").split(":")
                self.ld_template_list.append([key_snp, corresponding_snps.strip().split("\t")])
        print(F"\r[ {Fore.GREEN + 'x' + Style.RESET_ALL} ] Successfully read LD file (file: {src})    ")
    def createSpeciesMap(self, src):
        sys.stdout.write(F"\r[   ] Creating species map (file: {src})")
        sys.stdout.flush()
        # Open chromosome mapping file
        chromosome_mapping = open(src)
        
        for line in chromosome_mapping:
            # For each line, store (chr name, starting location)
            [chromosome, start] = line.replace("\n", "").split("\t")
            # Parse to int
            start = int(start)
            self.chromosome_structure.append((chromosome, start))
        
        # Potentially do: chromosome strucuture value
        self.chromosome_structure_val.append(0)
        for k in range(len(self.chromosome_structure)):
            self.chromosome_structure_val.append(self.chromosome_structure_val[-1] + self.chromosome_structure[k][1])
        # for (auto k = 0; k < RELI::chromosome_strucuture.size(); ++k){
            # chromosome_strucuture_val.push_back(chromosome_strucuture_val.back() + chromosome_strucuture.at(k).second);
            # //cout << chromosome_strucuture_val.back() << endl;
        # }
        print(F"\r[ {Fore.GREEN + 'x' + Style.RESET_ALL} ] Species map created (file: {src})  ")
    def readChIPSeqIndex(self, src, target):
        sys.stdout.write(F"\r[   ] Reading ChIP-seq index (file: {src})")
        sys.stdout.flush()
        index = 0
        # Open file, then for each line, read data info
        for line in open(src):
            if len(line) == 0:
                continue
            [datalabel, source, cell, tf, cell_label, pmid, group, ebv_status, species] = line.replace("\n", "").split("\t")
            self.chip_seq_index.append({
                "datalabel": datalabel,
                "source": source,
                "cell": cell,
                "tf": tf,
                "cell_label": cell_label,
                "pmid": pmid,
                "group": group,
                "ebv_status": ebv_status,
                "species": species
            })
            
            # While reading, check if this is the target
            if datalabel == target:
                # This is the target, find + read its ChIP-seq content
                self.target_data_location = os.path.join(self.directory, "ChIP-seq", target)
                self.target_data_index = index
            index += 1
        
        if self.target_data_index == -1:
            # Unable to locate target data
            print(F"\r[ {Fore.RED + 'x' + Style.RESET_ALL} ] Failure finding target data in ChIP-seq of ID {target}!")
            quitting(LINE())
        print(F"\r[ {Fore.GREEN + 'x' + Style.RESET_ALL} ] Parsed ChIP-seq index file (file: {src})")
    
    def loadNullData(self, src):
        # if (!RELI::snp_matching){	
            # in.ignore(bufferSize, '\n'); 
        # }
        sys.stdout.write(F"\r[   ] Reading null data (file: {src})")
        sys.stdout.flush()
        i = 0
        for line in open(src):
            # if (RELI::snp_matching){ 
                # RELI::binned_null_model_data.bin_map[atoi(linehandler(buffer).at(1).c_str())]->push_back(atoi(linehandler(buffer).at(0).c_str()));
            # }
            # else{
            self.null_model_data.append(int(line.replace("\n", "").split("\t")[0]))
            if i % 2000000 == 0:
                sys.stdout.write(F"\r[ \ ] Reading null data (file: {src})   ")
            elif i % 1000000 == 0:
                sys.stdout.write(F"\r[ - ] Reading null data (file: {src})   ")
            elif i % 500000 == 0:
                sys.stdout.write(F"\r[ / ] Reading null data (file: {src})   ")
            i += 1
            # }
        print(F"\r[ {Fore.GREEN + 'x' + Style.RESET_ALL} ] Parsed null data (file: {src})    ")
    
    def loadSNPFile(self, src):
        sys.stdout.write(F"\r[   ] Reading phenotype SNP data (file: {src})")
        sys.stdout.flush()
        
        for line in open(src):
            [snp_chr, snp_start, snp_end, snp_name] = line.replace("\n", "").split("\t")
            newsnp = {
                "snp_chr": snp_chr,
                "snp_start": int(snp_start),
                "snp_end": int(snp_end),
                "snp_name": snp_name
            }
            newsnp["length"] = newsnp["snp_end"] - newsnp["snp_start"]

            self.snp_vec.append(newsnp)

        print(F"\r[ {Fore.GREEN + 'x' + Style.RESET_ALL} ] Parsed phenotype SNP data (file: {src})    ")

    def loadDbSNPTable(self, src):
        sys.stdout.write(F"\r[   ] Reading dbSNP table (file: {src})")
        sys.stdout.flush()
        i = 0
        for line in open(src):
            linevec = line.replace("\n", "").split("\t")
            t = {
                "chr": linevec[0],
                "start": linevec[1],
                "end": linevec[2],
                "rsid": linevec[3],
                "obs_strand": linevec[4],
                "ref_allele": linevec[5],
                "alt_alleles": linevec[6],
                "type": linevec[7],
                "alt_allele_info": linevec[8],
                "alt_allele_freq": linevec[9]
            }

            self.snp_table_map[t["rsid"]] = t
            if i % 200000 == 0:
                sys.stdout.write(F"\r[ \ ] Reading dbSNP table (file: {src})")
            elif i % 100000 == 0:
                sys.stdout.write(F"\r[ - ] Reading dbSNP table (file: {src})")
            elif i % 50000 == 0:
                sys.stdout.write(F"\r[ / ] Reading dbSNP table (file: {src})")
            i += 1
            
        print(F"\r[ {Fore.GREEN + 'x' + Style.RESET_ALL} ] Parsed dbSNP table (file: {src})    ")

    def loadLDSNPs(self, src):
        sys.stdout.write(F"\r[   ] Reading LD SNP file (file: {src})")
        sys.stdout.flush()
        if src != None:
            # Read through the LD file
            self.read_ld(src)
            for entry in self.ld_template_list:
                # Entry is in format [keySNP, mySNP]
                newld = {
                    # Find matching keySNP in SNP_vec_temp
                    "keySNP": next(x for x in self.SNP_vec_temp if x["snp_name"] == entry[0]),
                    "mySNP": [],
                    "dis2keySNP": []
                }
                # Iterate through mySNP entries
                for i in range(len(entry[1])):
                    # Find the SNP entry with this value at this location i
                    for (ii, snpit) in enumerate(self.SNP_vec_temp):
                        if entry[1][i] == snpit["snp_name"]:
                            # We found it! Push into mySNP + erase from tmp vector
                            newld["mySNP"].append(snpit)
                            del self.SNP_vec_temp[ii]
                            break
                self.ld_list.append(newld)
        # Go through all leftover values + push into LD list
        for leftover_entry in self.SNP_vec_temp:
            self.ld_list.append({
                "keySNP": leftover_entry,
                "mySNP": [leftover_entry],
                "dis2keySNP": []
            })
        
        for ldit in self.ld_list:
            # Iterate through mySNP
            for snpit in ldit["mySNP"]:
                ldit["dis2keySNP"].append(snpit["snp_end"] - ldit["keySNP"]["snp_end"])
        # print(self.ld_list)
        print(F"\r[ {Fore.GREEN + 'x' + Style.RESET_ALL} ] Parsed LD SNP table (file: {src})    ")

    def SNPfit(self, LD_A, tempSNP_A, int_A, inVec, 
	inVec2):
        pos_set = []
        neg_set = []
        max_diff = 0
        for it in LD_A["dis2keySNP"]:
            if it > 0:
                pos_set.append(it)
            else:
                neg_set.append(it)
        if len(pos_set) > 0:
            max_diff = max(pos_set)
        if len(neg_set) > 0:
            max_diff = max(max_diff, abs(min(neg_set)))
        
        okay = False
        for k in range(len(inVec)):
            if (int_A - max_diff - tempSNP_A["length"] >= inVec2[k] and int_A + max_diff + tempSNP_A["length"] <= inVec2[k + 1]):
                tempSNP_A["snp_chr"] = inVec[k][0]
                tempSNP_A["snp_end"] = int_A + math.floor(tempSNP_A["length"] / 2) - inVec2[k]
                tempSNP_A["snp_start"] = tempSNP_A["snp_end"] - tempSNP_A["length"]
                okay = True
                break
        return okay
    
    def overlapping(self, SNPvecA, bedvecA, rsid_collector, _iter_number):
        temp_snp_vec = SNPvecA
        in_LD_unique_key_collector = []
        k = 0
        temp_snp_vec.sort(key=lambda x: (x["snp_chr"],x["snp_start"]))
        prev_chr = "chr0"
        
        for snpit in temp_snp_vec:
            if snpit["snp_chr"] == prev_chr:
                t = max(max(k - self.lookback, 0), self.target_bed_index[snpit["snp_chr"]])
                k = t
                while k < len(bedvecA):
                    if bedvecA[k]["bed_chr"] == snpit["snp_chr"] and overlap_exists(bedvecA[k], snpit):
                        in_LD_unique_key_collector.append(snpit["inherited_unique_key_from_LD"])
                        if (_iter_number == 0):
                            rsid_collector.append(snpit["snp_name"])
                        break
                    if (bedvecA[k]["bed_start"] >= snpit["snp_end"] or bedvecA[k]["bed_chr"] != snpit["snp_chr"]):
                        break
                    k += 1
            else:
                prev_chr = snpit["snp_chr"]
                k = self.target_bed_index[snpit["snp_chr"]]
                while k < len(bedvecA):
                    if bedvecA[k]["bed_chr"] == snpit["snp_chr"] and overlap_exists(bedvecA[k], snpit):
                        in_LD_unique_key_collector.append(snpit["inherited_unique_key_from_LD"])
                        if (_iter_number == 0):
                            rsid_collector.append(snpit["snp_name"])
                        break
                    
                    if bedvecA[k]["bed_start"] >= snpit["snp_end"] or bedvecA[k]["bed_chr"] != snpit["snp_chr"]:
                        break
                    k += 1
        
        return in_LD_unique_key_collector
    
    def sim(self):
        for i in range(self.num_reps + 1):
            self.ld_sim_vec = []
            start_pt = self.ld_list[0]
            if i == 0:
                # First time around
                for (ind, ldit) in enumerate(self.ld_list):
                    ld_sim = {
                        "unique_key": ind,
                        "mySNP": ldit["mySNP"]
                    }
                    # Iter through SNP list
                    for tSNP in ld_sim["mySNP"]:
                        tSNP["inherited_unique_key_from_LD"] = ld_sim["unique_key"]
                    
                    ld_sim["dis2keySNP"] = ldit["dis2keySNP"]
                    ld_sim["keySNP"] = ldit["keySNP"]
                    ld_sim["overlap_sim"] = False
                    self.ld_sim_vec.append(ld_sim)
            else:
                for (ind, ldit) in enumerate(self.ld_list):
                    ld_sim = { "unique_key": ind, "mySNP": [] }
                    t_index = 0
                    t_key_snp = { "length": ldit["keySNP"]["length"] }
                    datagood = False
                    # Add if SNP matching
                    """
                    if (RELI::snp_matching){
                        std::uniform_int_distribution<unsigned int> distGen(0, (RELI::binned_null_model_data.bin_map[LDit->keySNP._MAF_Bin]->size() - 1));
                        while (datagood != true){
                            tIndex = distGen(randSeed);
                            datagood = RELI::SNPfit(*LDit, 
                                tKeySNP, 
                                RELI::binned_null_model_data.bin_map[LDit->keySNP._MAF_Bin]->at(tIndex),
                                RELI::chromosome_strucuture, 
                                RELI::chromosome_strucuture_val, 
                                1);	
                        }
                    }
                    """
                    # random.randint(0, len(self.null_model_data) - 1)
                    while not datagood:
                        t_index = random.randint(0, len(self.null_model_data) - 1)
                        datagood = self.SNPfit(ldit,
                            t_key_snp,
                            self.null_model_data[t_index],
                            self.chromosome_structure,
                            self.chromosome_structure_val
                        )
                    for dit in ldit["dis2keySNP"]:
                        t_snp = {}
                        # Fill t_snp with information from the key SNP
                        t_snp["length"] = t_key_snp["length"]
                        t_snp["snp_chr"] = t_key_snp["snp_chr"]
                        t_snp["snp_end"] = t_key_snp["snp_end"] + dit
                        t_snp["snp_start"] = t_snp["snp_end"] - t_snp["length"]
                        
                        t_snp["inherited_unique_key_from_LD"] = ld_sim["unique_key"]
                        ld_sim["mySNP"].append(t_snp)
                    
                    ld_sim["overlap_sim"] = False
                    self.ld_sim_vec.append(ld_sim)
            self.SNP_vec_temp = []
            for ldsimit in self.ld_sim_vec:
                for snpit in ldsimit["mySNP"]:
                    self.SNP_vec_temp.append(snpit)
            ld_key_collector = self.overlapping(self.SNP_vec_temp, self.target_bed_vec, self.overlapped_rsids, i)
            ld_key_collector.sort()
            # Now count all unique values
            ld_unique_key_collector = set(ld_key_collector)
            if i % 500 == 0:
                print(F"{float(i)/float(self.num_reps)*100}% finished.")
            print(F"Current iteration: {i}, current intersection: {len(ld_unique_key_collector)}")
            self.stats_vec.append(float(len(ld_unique_key_collector)))

    def cal_stats(self, model_mode):
        sys.stdout.write(F"\r[   ] Finished analysis, parsing statistics (MODE: {model_mode})...")
        sys.stdout.flush()
        # Mean value
        self.mu = sum(self.stats_vec) / float(len(self.stats_vec))
        
        temp = 0
        for it in self.stats_vec:
            temp += (it - self.mu)*(it - self.mu)

        # Sample standard deviation
        self.sd = math.sqrt(temp / float(len(self.stats_vec) - 1))
        if self.sd == 0 or self.stats_vec[0] < math.ceil(float(len(self.ld_list))*self.sig_pct):
            self.zscore = 0
        else:
            self.zscore = (self.stats_vec[0] - self.mu) / self.sd
        
        if model_mode is NORMAL:
            # Using normal distribution
            self.pval = scipy.stats.norm.sf(self.zscore)
            self.corr_pval = min(self.pval*self.corr_muliplier, 1.0)
        elif model_mode is EMPIRICAL:
            real_obs = self.stats_vec[0]
            tvec = self.stats_vec
            tvec.sort()
            
            # Number of instances >= the first value
            greater_or_equal_instance = len(list(filter(tvec, lambda x: x >= real_obs)))
            self.pval = scipy.stats.norm.sf(self.zscore)
            self.corr_pval = float(greater_or_equal_instance) / float(len(self.stats_vec))
        elif model_mode is DEFAULT:
            self.pval = scipy.stats.norm.sf(self.zscore)
            self.corr_pval = min(self.pval*self.corr_muliplier, 1.0)
        print(F"\r[ {Fore.GREEN + 'x' + Style.RESET_ALL} ] Successfully parsed statistics. All done!          ")
    def interpret_output(self, model):
        # Calculate output results
        self.cal_stats(model)
        # Create a folder, if it doesn't already exists
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        print_output = open(os.path.join(self.output_dir, "summary.txt"), 'w+')
        print_output.write(F"{self.phenotype_name} Analysis:\n")
        print_output.write(F"\tTF: {self.chip_seq_index[self.target_data_index]['tf']}\n")
        print_output.write(F"\tCell Label: {self.chip_seq_index[self.target_data_index]['cell_label']}\n")
        print_output.write(F"\tCell: {self.chip_seq_index[self.target_data_index]['cell']}\n")
        print_output.write(F"\tIntersect: {self.stats_vec[0]}\n")
        print_output.write(F"\tTotal: {len(self.ld_list)}\n")
        print_output.write(F"\tRatio: {self.stats_vec[0] / len(self.ld_list)}\n")
        print_output.write(F"\tMean: {self.mu}\n")
        print_output.write(F"\tStd: {self.sd}\n")
        print_output.write(F"\tZ-Score: {self.zscore}\n")
        print_output.write(F"\tRelative Risk: {self.stats_vec[0] / self.mu if self.mu != 0 else 0}\n")
        print_output.write(F"\tP Value: {self.pval}\n")
        print_output.write(F"\tCorrected P Value: {self.corr_pval}")
        
        print_output.close()
        
        rsid_list = open(os.path.join(self.output_dir, "result.overlaps"), "w+")
        for overlaps in self.stats_vec:
            rsid_list.write(str(int(overlaps)) + "\n")
        rsid_list.close()
        overlapped = list(set(self.overlapped_rsids))
        f = open(os.path.join(self.output_dir, "result.rsids"), "w+")
        f.write("\n".join(overlapped))
        f.close()

    def __init__(
        self,
        phenotype_name,
        snp_file,
        num_reps,
        ld_file,
        chipseq_index,
        target,
        directory = "sample_data",
        null = "sample_data/Null/CommonSNP_MAFmatch",
        dbsnp_index = "sample_data/SNPtable/SNPtable",
        given_species = "hg19.txt",
        output_dir = "output"
    ):
        # Seed this
        seed = int(time.time())
        print(F"Starting RELI:\n\tSeed: {seed}\n\tSNP File: {snp_file}")
        # Read and initialize all the information given
        self.num_reps = num_reps
        self.ld_list = []
        self.phenotype_name = phenotype_name
        # Information about each ChIPseq entry given
        self.chip_seq_index = []
        # Index where the target data is stored in chip_seq_index
        self.target_data_index = -1
        # Tuple of all chromosome starting locations
        self.chromosome_structure = []
        # 
        self.chromosome_structure_val = []
        # Directory where default information is stored
        self.directory = directory
        # Where null model information is stored
        self.null_model_data = []
        # List of phenotype SNP entries
        self.snp_vec = []
        # List of dbSNP entries, indexed by RSID
        self.snp_table_map = {}
        # Used by read_ld to store intermediary informtion
        self.ld_template_list = []
        self.stats_vec = []
        self.sig_pct = 0.05
        self.corr_muliplier = 1
        # TODO: understand this
        self.lookback = 50
        # List of overlapped RSIDs found
        self.overlapped_rsids = []
        
        # Load genome structure - by default human species
        self.createSpeciesMap(given_species)
        # Read ChIP-seq index file and set target ChIP-seq file
        self.readChIPSeqIndex(chipseq_index, target)
        # Init target ChIP-seq file object
        # RELI::target_bed_file TBF;
        
        # Load target ChIP-seq file (target_data_location set by readChIPSeqIndex)
        tbf = BEDSample(self.target_data_location)
        # Store target BED vector + index in structure
        self.target_bed_vec = tbf.data
        # Contains the location of where each chromosome starts
        self.target_bed_index = tbf.index
        # Load null model
        # (take list of common SNPs)
        self.loadNullData(null)
        # Load phenotype SNP file
        self.loadSNPFile(snp_file)
        # Load dbSNP table
        self.loadDbSNPTable(dbsnp_index)
        # Extract snp bin info from dbsnp table - doesn't do anything if snp_matching disabled
        # self.extractSNPInfo()
        # RELIinstance->extract_snp_info(RELIinstance->ATGCmap);
        # Copy over SNP data for loading the LD structure
        self.SNP_vec_temp = self.snp_vec
        # Load phenotype LD structure
        self.loadLDSNPs(ld_file)
    
        # Begin the simulation!
        self.sim()
        
        # Display output results
        self.output_dir = output_dir
        self.interpret_output(NORMAL)

instance = RELI(
    "SLE",
    "example/SLE_EU.snp",
    2000,
    "example/SLE_EU.ld",
    "sample_data/ChIPseq.index",
    "hg19_0302",
    given_species = "sample_data/GenomeBuild/hg19.txt"
)
