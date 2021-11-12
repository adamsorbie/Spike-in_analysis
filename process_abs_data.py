#!/usr/bin/env python
from Bio import SeqIO
from sys import platform
import argparse
import pandas as pd
import numpy as np
import subprocess
import os


## to-do 
# include on linux flag which will run more sensitive bowtie2 alignment 
# remove ng spike requirement
# double check against core facility script 
# add more error handling and warnings 
# remove irrelevant documentation 



## helper functions
# normalise otu table by minimum sum scaling

	
def normalise(otu, rel=False):
    """

    :param  pandas.Dataframe otu: OTU table to be normalised
    :param rel: Whether or not to return OTU normalised by relative abundance
    :return: OTU table normalised by minimum sum scaling or relative abundance
    """
    if otu.shape[1] == otu.select_dtypes(include=np.number).shape[1]:
        sum_row = otu.sum()
        if rel:
            rel_abund = otu * 100 / sum_row
            return rel_abund
        else:
            min_sum = min(sum_row)
            normalised = otu / sum_row * min_sum
            return normalised
    else:
        print("Non-numeric values are not allowed")


# make directories from list
def make_dir(dirnames):
    """
    :param list dirnames: list of directories to create
    :return:
    """
    for i in dirnames:
        try:
            os.mkdir(i)
        except FileExistsError:
            print("Directory ", i, " already exists")


# store otu table and taxonomy column as separate variable
def drop_taxonomy(otu, keep_col_as_var=True):
    """

    :param pandas.Dataframe otu: OTU table
    :param keep_col_as_var: Store taxonomy column as a variable
    :return: OTU table and taxonomy column as pandas.Series
    """
    if "taxonomy" in otu.columns:
        if keep_col_as_var:
            tax_col = otu[["taxonomy"]]
            otus_only = otu.drop("taxonomy", axis=1)
            return otus_only, tax_col
        else:
            otus_only = otu.drop("taxonomy", axis=1)
            return otus_only
    else:
        print("taxonomy column does not exist")
        return otu, None

    
def check_order(otu, mapping):
    """

    :param  pandas.Dataframe otu: OTU table
    :param pandas.Dataframe mapping: mapping file
    :return: OTU table sorted by column name and mapping sorted by index
    """
    # drop taxonomy, store as variable and remove water columns
    otu_notax, tax_col =  remove_water_controls(otu, mapping, norm=False, return_taxa_as_var=True)
    # transpose mapping file and order by columns
    mapping_t = mapping.T
    # order mapping file by OTU cols, transpose and re-add taxa col 
    mapping_ordered = mapping_t[otu_notax.columns] 
    mapping_out = mapping_ordered.T
    otu_out = pd.concat([otu_notax, tax_col], axis=1)
    return otu_out, mapping_out


def generate_analysis_record(list_np_array, index_df, col_names):
    """

    :param list of numpy.ndarray list_np_array: list of numpy arrays to convert into columns
    :param pandas.Dataframe.index index_df: index object to be used as dataframe index
    :param list col_names: list of column names for the output dataframe
    :return: pandas.Dataframe, with each numpy array as a column
    """
    # hacky way to generate sorted index
    sorted_df = index_df.sort_index()
    sorted_index = sorted_df.index
    # make dataframe from list of vectors
    record_df = pd.DataFrame(list_np_array)
    # transpose to vectors are in columns
    record_df_out = record_df.T

    # set col names and add sorted index
    record_df_out.columns = col_names
    record_df_out.index = sorted_index
    return record_df_out


def remove_water_controls(otu, mapping, norm=True, return_taxa_as_var=False, drop_taxa_col=False):
    """


    :param pandas.DataFrame otu: OTU table
    :param pandas.DataFrame mapping: mapping file
    :param bool norm: Whether to return normalised OTU table or not
    :param bool return_taxa_as_var: Whether to return taxonomy column as variable
    :param bool drop_taxa_col: Drop taxonomy column completely if true
    :return: OTU table with water and other controls removed
    """
    otu, tax_col = drop_taxonomy(otu)
    # should add some check if there is at least some overlap in mapping index and OTU columns here 
    mapping_samples = mapping.index.tolist()
    water_cols = [col for col in otu.columns if col not in mapping_samples]
    otu_out = otu.drop(water_cols, axis=1)
    # this is a little messy
    if norm:
        otu_out_norm = normalise(otu_out)
        if return_taxa_as_var:
            return otu_out_norm, tax_col
        else:
            if drop_taxa_col:
                return otu_out_norm
            else:
                otu_out_norm["taxonomy"] = tax_col
                return otu_out_norm
    else:
        if return_taxa_as_var:
            return otu_out, tax_col
        else:
            if drop_taxa_col:
                return otu_out
            else:
                otu_out["taxonomy"] = tax_col
                return otu_out


## functions for finding potential spikes
# find potential spikes in fasta file based on taxonomy (spikes should only be classified as bacteria)
def find_spikes(otu, fasta):
    """

    find otus with no taxonomy (suggested spikes)
    and filter fasta file for these
    :param pandas.Dataframe otu: OTU table
    :param str fasta: fasta filename
    :return: Seq object containing putative spikes
    """
    if "taxonomy" in otu.columns:
        if not otu[otu.taxonomy == "Bacteria;;;;;;"].empty:
            otu_put_spikes = otu[otu.taxonomy == "Bacteria;;;;;;"]
            spike_otus = otu_put_spikes.index.tolist()
            otu_seqs = list(r for r in SeqIO.parse(fasta, "fasta") if r.id in spike_otus)
            SeqIO.write(otu_seqs, "putative_spikes.fasta", "fasta")
        else:
            print("'Bacteria;;;;;' not found, does your input OTU contain spikes?'")
    else:
        print("No taxonomy column found")


# function to call blast from command line
def blast_seqs(query, subject, out_file):
    """ Runs BLAST, searching for spike sequences in OTU fasta file.

    :param str query: Query fasta filename
    :param str subject: Subject fasta filename
    :param str out_file: Name of output text file containing hits
    :return: None
    """
    blast_cmd = "blastn -task megablast -query " + str(query) + "  -subject " + str(subject) + " -outfmt 6 -out " + str(
        out_file) + " "
    subprocess.call(blast_cmd, shell=True)
    
# def run_bowtie2(query, subject, out_file):
#     """ Runs bowtie2, searching for spike sequences in OTU fasta file (Linux only).

#     :param str query: Query fasta filename
#     :param str subject: Subject fasta filename
#     :param str out_file: Name of output text file containing hits
#     :return: None
#     """
#     if platform not in ["linux", "darwin"]:
#         print("Bowtie2 is only supported on Linux and Mac")
#     else: 
#         bowtie2_cmd = ["bowtie2", "-x", spikes_ref_name, "-1", fastq_R1, "-2", fastq_R2, "--al-conc", fastq_aligned, "--un-conc",
#            fastq_unaligned]


# filter blast results for hits above 98% identity
def filter_blast_results(blast_out):
    """

    :param str blast_out: filename of BLAST results
    :return: OTUs > 98% identity with Spikes and eval < 1e-50
    """
    blast_res = pd.read_csv(blast_out, sep="\t", header=None)
    if not blast_res[blast_res.iloc[:, 2] > 98].empty:
        spike_hits = blast_res[blast_res.iloc[:, 2] > 98]
    else:
        print("No hits with significant identity found")
        return None
    # ideally this goes in blastn call but can't figure out how to properly format as float
    if not spike_hits[spike_hits.iloc[:, 10] < 1e-50].empty:
        spike_hits_conf = spike_hits[spike_hits.iloc[:, 10] < 1e-50]
        spike_hits_conf.to_csv("confirmed-spikes.csv")
    else:
        print("No hits below evalue threshold of 1e-50 found")
        return None
    return spike_hits_conf[0].unique().tolist()


# functions to remove spikes
def rem_spikes(fasta, list_spikes):
    """

    remove spikes from input fasta file
    :param str fasta: Filename of fasta containing OTU sequences
    :param str list_spikes: list of confirmed spike OTUs (output of filter_blast_results)
    :return: Write filtered Sequences to fasta file
    """
    if list_spikes:
        otu_seqs = list(r for r in SeqIO.parse(fasta, "fasta") if r.id not in list_spikes)
    else:
        print("Spikes not in input file")
        return None
    SeqIO.write(otu_seqs, "OTUs-Seqs-filt.fasta", "fasta")


# functions to format and output relative abundance tables
def generate_rel_otu(otu, mapping, list_spikes):
    """
    generate relative abundance OTU tables for analysis
    :param pandas.Dataframe otu: input OTU table
    :param pandas.Dataframe mapping: mapping file with only samples for analysis
    :param list list_spikes: list of confirmed spike OTUs
    :return: filtered OTU tables (un-normalised and normalised)
    """
    # need to think more about error handling here
    # remove spike otus from input OTU
    otu_no_spikes = otu[~otu.index.isin(list_spikes)]
    # remove water controls from above returning un-normalised OTU and taxonomy col as df
    otu_samples_only, tax_col = remove_water_controls(otu_no_spikes, mapping, norm=False, return_taxa_as_var=True)
    
	# normalise output of above call and add taxonomy col 
    otu_norm = normalise(otu_samples_only)
    otu_norm["taxonomy"] = tax_col
    # add taxonomy call to un-normalised OTU
    otu_samples_only["taxonomy"] = tax_col
    return otu_samples_only, otu_norm

### edit generate and calculate abs abund to include df (multiply spike mean by this)
# functions to calculate absolute abundance and write out
def calc_abs_abund(otu, list_spikes, mapping, weight_col, df=None):
    """

    :param pandas.Dataframe otu: OTU table
    :param list list_spikes: list of spike OTUs
    :param pandas.Dataframe mapping: mapping file
    :param str weight_col: column in mapping file with sample weights
    :return: numpy array of bacterial reads per mg
    """
	# add check which ensures otu and mapping are in same order 
    # get total spike reads per sample
    spike_rps = np.array(otu[otu.index.isin(list_spikes)].sum())
    
    if df is not None:
        spike_mean = np.sum(spike_rps) / len(spike_rps) * df
    else:
        spike_mean = np.sum(spike_rps) / len(spike_rps)
    # get total bacterial reads per sample for records
    bacterial_rps = np.array(otu[~otu.index.isin(list_spikes)].sum())
    # remove spikes from otu table
    otu_spikes_rem = otu.drop(list_spikes)
    # store sample weights as vector
    sample_weights = mapping[weight_col].to_numpy()
    # calculate abs abund
    abs_otu = (otu_spikes_rem * spike_mean) / (spike_rps * sample_weights)
    # list of calculations to keep a record of
    list_np = [spike_rps, bacterial_rps, sample_weights]
    # column names for output data file
    col_list = ["Spike_RPS", "Bacterial_RPS", "Sample_weights"]
    # use generate_analysis_record function to create dataframe
    calculations_record = generate_analysis_record(list_np, abs_otu.T, col_list)

    return abs_otu, calculations_record


def generate_abs_abund(original_otu, list_spikes, tax_col, **kwargs):
    """

    :param pandas.Dataframe original_otu: un-normalised OTU table
    :param list_spikes: list of spike OTUs
    :param pandas.Dataframe tax_col: taxonomy column
    :param `**kwargs`: Keyword arguments are passed to calc_abs_abund
    :return: Absolute abundance OTU table 
    """
    # full calculation: Bacterial RPS / (Sum spike reads/adjust to * spike/spike)/sample weight (current more readable)
    abs_otu, calc_records = calc_abs_abund(original_otu, list_spikes, **kwargs)
    
    tax_col_spikes_rem = tax_col.drop(list_spikes)


    abs_otu["taxonomy"] = tax_col_spikes_rem
    return abs_otu, calc_records


## main function

def main():
    """
    parse cli arguments and write output to file
    :return:
    """
    parser = argparse.ArgumentParser(description="DESCRIPTION\n"
                                    "This script allows for automatic processing of 16S sequencing data \n"
                                    "which has been 'spiked' with a known amount of synthetic sequences \n"
                                    "to allow for the calculation of absolute abundance.\n"
                                    "\n\n==========================BASIC USAGE==========================\n"
                                    "\n$ process_abs_data.py -f in.fna -otu otu.tab -m meta.tab\n"
                                    "-s spikes.fna", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-f", "--fasta", required=True, type=str, help="name of OTU fasta file")
    parser.add_argument("-otu", "--otu_table", required=True, type=str, help="original unaltered otu")
    parser.add_argument("-m", "--mapping", required=True, type=str, help="mapping file")
    parser.add_argument("-s", "--spikes", required=True, type=str, help="fasta file containing spike sequences")
    parser.add_argument("-w", "--weight_col", required=True, type=str,
                        help="column in mapping file with sample weights")
    parser.add_argument("-df", "--dilution_factor", required=False, default="", type=int, help="dilution factor of spike")
    args = parser.parse_args()

    # load required files
    otu = pd.read_csv(args.otu_table, sep="\t", header=0, index_col=0)
    mapping = pd.read_csv(args.mapping, sep="\t", header=0, index_col=0)
	
	# sort by sample name
    otu, mapping = check_order(otu, mapping)

    # make directories
    dirs = ["relative_abund", "absolute_abund", "log"]
    make_dir(dirs)

    # find spikes
    find_spikes(otu, args.fasta)
    blast_seqs(query="putative_spikes.fasta", subject=args.spikes, out_file="results.txt")
    list_spikes = filter_blast_results("results.txt")

    # remove confirmed spikes from fasta file
    rem_spikes(args.fasta, list_spikes)
    os.system("cp mapping_file.tab OTUs-Seqs-filt.fasta relative_abund/")
    os.system("cp mapping_file.tab OTUs-Seqs-filt.fasta absolute_abund/")

    # generate relative abundance OTU table -
    rel_otu, rel_otu_norm = generate_rel_otu(otu, mapping, list_spikes)
    rel_otu.to_csv("relative_abund/OTUs-Table-Filt.tab", sep="\t")
    rel_otu_norm.to_csv("relative_abund/OTUs-Table-Filt-Norm.tab", sep="\t")

    # create normalised OTU with spikes for calculation only
    norm_spikes_otu, tax_col = remove_water_controls(otu, mapping, norm=True, return_taxa_as_var=True)
    # original otu for spike removal, post-hoc normalisation and abs abund calc
    original_otu = remove_water_controls(otu, mapping, norm=False, return_taxa_as_var=False, drop_taxa_col=True)
    # absolute abundance
    abs_abundance, calc_records = generate_abs_abund(original_otu,
                                                    list_spikes,
                                                    tax_col,
                                                    mapping=mapping,
                                                    weight_col=args.weight_col,
                                                    df=args.dilution_factor)

    # named args passed to calc
    abs_abundance.to_csv("absolute_abund/OTUs-Table-Abs-Norm.tab", sep="\t")

    # write out calculations
    calc_records.to_csv("log/calculations_record.csv")

    # tidy up folder
    os.system("mv putative_spikes.fasta results.txt confirmed-spikes.csv log/")
    os.system("rm -f OTUs-Seqs-filt.fasta")


if __name__ == '__main__':
    main()
