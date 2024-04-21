<h1>RELI</h1>

INPUT
------

1,545 SNPs in input index file
14,063,450 SNPs given

90 RSIDs, 2000 overlaps

Constants involved:

- public_ver_snp_fname - name of the phenotype SNP file being used.
- ldfile - name of the phenotype LD structure
- snp_matching - whether or not SNP matching is enabled or disabled.
- ldfile_flag - if LD file was specified (SLE_EU.ld)

Steps (RELI.cpp main with references from RELI_impl):

1. Load genome structure w/ default (hg19 species)
2. Create pair<chromosome #, location # in genome> + push into chromosome_structure vector if default, else push from file
3. Read ChIP-seq index file (public_ver_data_index_fname)
    - Read line by line and push each element into dataindexvec (store datalabel, source, cell, tf, cell_label, pmid, group, ebv_status, and species)
4. Set target ChIP-seq file
    - Search through all of the ChIP-seq files for the line corresponding to the public_ver_target_label. This is by default hg19_0302 (a human person that tested EBV+; this is an analysis of the binding sites of the EBNA2 protein to this person’s genome). This is put into public_ver_target_data_fname.
5. Init target file object, load target ChIP-seq file
    - Initialize file object of type RELI::target_bed_file
    - Take the target ChIP-seq file; for each line, set bed_chr to the chromosome, then take the start and ends and put into bed_start and bed_end, set length to bed_end - bed_start and push length lengthvec and the bed3col into myData.
    - create target ChIP-seq file index - start at the first line’s bed_chr, and set that string, prev_chr, at the index to 0. Then go along the data and if the line isn’t the same as prev_chr, set it to prev_chr and take the distance between that point and the beginning of the data.
6. Set the targetbedinfilevec to myData (the target ChIP-seq file) and the targetbedinfilevecindex_start to index.
7. Get null model data from CommonSNP_MAF (Major Allele Frequency) match. This contains a list of all common SNPs (possibly); it’s read by setting the nullmodelinfilename to the name of the file, then reading the file line-by-line.
    - If matching is enabled, take the number of SNPs matched (?), access bin_map at that index, and push the location of the SNPs into bin_map at this index.
    - If matching is not enabled (this is the default), push all locations into one bin, bin0.
8. Load phenotype SNP file - go through each line of the public_ver_snp_fname file and for each store the name of the chromosome (snp_chr), the start, end, and name, and push into SNP_vec
9. Load DB SNP table - go through each line of the public_ver_snp_table_fname and store the chromosome name (chr), start, end, rsid, obs_strand (?, either +/-), ref_allele, alt_alleles, type, alt_allele_info, alt_allele_freq; store into the snptablemap, indexed by rsid
10. Extract SNP info - If SNP matching, iterate through SNP_vec (containing the phenotype SNP file), then see if the SNP name is in the snptablemap (dbSNP table). If it is, update the SNP_vec with the obs_strand,_ref_allele, and snp_type (vs type) *TODO*
11. Copy over the SNP data into SNP_vec_temp for loading the LD structure.
12. Load phenotype LD structure file (load_ld_snps).
    - Read LD file: for each line, push LD template, and put each rsid, ex. rs13023380, into mySNP, and keySNP as the key in LD_template_vec.
    - Iterate over LD_template_vec and, for each entry, search through SNP_vec_temp to find keySNP + place into new instance. **
    -
13. Perform the actual simulation repmax number of times:
    a. Iterate over LD_vec


