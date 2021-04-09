### goal of this script is to simulate an approach to estimate proteomic profiles from metaproteomic data

### notation:
### p = unique peptide
### k = coarse grained group
### j = organism
### r = redundant peptide

import numpy as np
import math
from random import randrange
import random
from Bio.Alphabet import IUPAC
from pyteomics import mass
from pandas import DataFrame
import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns

##### GENERATIVE MODEL

## within each j, draw from k different groups which have more or less sequence divesity
# make k sequence banks of different sizes
def generate_k_seq_banks(number_of_seqs):
    '''
    Generate sequence banks of different compositions. Bigger sequence banks
    implies a more diverse coarse grained pool.
    '''
    amino_acids = IUPAC.IUPACProtein.letters.upper()
    sequence_bank = []
    # a coarse grained pool of size number_of_seqs
    for seq in range(number_of_seqs):
        # choose a length between 5 and 15
        seq_length = randrange(10) + 5
        # make a sequence of length
        my_seq = ''.join(random.choice(amino_acids) for i in range(seq_length))
        # append the sequence to a list
        sequence_bank.append(my_seq)
    return(sequence_bank)

def generate_all_seq_banks(number_of_seqs_vec):
    '''
    Take the function above and do that for many k values
    '''
    all_seq_banks = []
    for k in number_of_seqs_vec:
        seq_bank_k = generate_k_seq_banks(number_of_seqs = k)
        all_seq_banks.append(seq_bank_k)
    return(all_seq_banks)

## generate j proteomic profiles in the form of peptide sequences.
def generate_organism_profile(all_seq_banks_out, n_seqs):
    '''
    Take the all seq banks function output, and subsampled sequences from it
    '''
    peptide_profile = []
    protein_pool_id = []
    i = 0

    for k_group in all_seq_banks_out:
        sub_group = k_group[:]
        for n_seq in range(n_seqs):
            random_peptide = random.choice(sub_group)
            sub_group.remove(random_peptide)
            peptide_profile.append(random_peptide)
            protein_pool_id.append(i)
        i += 1
#     print(protein_pool_id)
    return([peptide_profile, protein_pool_id])

def generate_organism_profiles(all_seq_banks_out, n_seqs, n_orgs):
    '''
    generate many organism peptide profiles, returning a list of lists
    '''
    organism_profiles = []
    organism_protein_pool_ids = []

    for organism in range(n_orgs):
        organism_profile_i = generate_organism_profile(all_seq_banks_out = all_seq_banks_out,
                                                       n_seqs = n_seqs)
        organism_profiles.append(organism_profile_i[0])
        organism_protein_pool_ids.append(organism_profile_i[1])

    return([organism_profiles, organism_protein_pool_ids])

def generate_tax_abundance(n_orgs):
    '''
    generate taxon specific abundances which will then be used to multiply peptide abundances for a given taxon
    '''
    tax_abundances = []
    for organism in range(n_orgs):
        tax_abundances.append(np.random.gamma(1, 10))
    return(tax_abundances)

def generate_pool_abundance(n_pools):
    '''
    generate taxon specific abundances which will then be used to multiply peptide abundances for a given taxon
    '''
    pool_abundances = []
    for pool in range(n_pools):
        pool_abundances.append(np.random.gamma(1, 10))
    return(pool_abundances)

def generate_peptide_abundance(organism_profiles_out,
                               protein_pool_abundance,
                               tax_abundances_out):
    '''
    go through all the peptide abundances and draw values from a gamma distribution
    '''
    organism_peptide_abundances = []

    # go through every organism
    for organism_pept in range(len(organism_profiles_out)):
        peptide_profile_i = []
        organism_peptides_i = organism_profiles_out[organism_pept]
        # go through every peptide for a given organism
        for protein_pool_i in protein_pool_abundance:
            for peptide_i in range(len(organism_peptides_i)//len(protein_pool_abundance)):
                # generate a peptide abundance, but multiply by the abundance of a given organism/protein_pool
                 # for every protein pool there is a separate abundance factor
                peptide_profile_i.append(np.random.gamma(0.15, 10)*tax_abundances_out[organism_pept]*protein_pool_i)
        organism_peptide_abundances.append(peptide_profile_i)

    return(organism_peptide_abundances)

def write_dataset(n_orgs, n_seqs,
                  number_of_seqs_vec = [250, 1000, 2500, 10000],
                  rt_low = 0,
                  rt_upp = 120):

    all_seq_banks_out = generate_all_seq_banks(number_of_seqs_vec = number_of_seqs_vec)
    organism_profiles_out = generate_organism_profiles(all_seq_banks_out = all_seq_banks_out,
                                                       n_seqs = n_seqs,
                                                       n_orgs = n_orgs)
    tax_abundances_out = generate_tax_abundance(n_orgs)
    pool_abundances_out = generate_pool_abundance(n_pools = len(number_of_seqs_vec))

    peptide_abundances = generate_peptide_abundance(organism_profiles_out = organism_profiles_out[0],
                                                    protein_pool_abundance = pool_abundances_out,
                                                    tax_abundances_out = tax_abundances_out)

    flat_peptide_abundances = [item for sublist in peptide_abundances for item in sublist]
    flat_organism_profiles = [item for sublist in organism_profiles_out[0] for item in sublist]
    flat_organism_protein_pools = [item for sublist in organism_profiles_out[1] for item in sublist]
#     print(flat_organism_protein_pools)
    # make a list of organism IDs
    # number of entires in total
    num_entries = n_orgs*n_seqs*len(number_of_seqs_vec)
    organism_index = []
    for organism in range(n_orgs):
        organism_index.extend([organism]*n_seqs*len(number_of_seqs_vec))

    # generate an identifier for the protein pool
#     protein_pool_index = []
#     for protein_pool in range(len(number_of_seqs_vec)):
#         protein_pool_index.extend([protein_pool]*n_seqs*n_orgs)

    # generating truncated uniform retention times
    retention_times = np.random.uniform(rt_low, rt_upp, num_entries)

    # calculating peptide mass
    peptide_mass = []
    for peptide in flat_organism_profiles:
        peptide_mass_i = mass.calculate_mass(sequence = peptide)
        peptide_mass.append(peptide_mass_i/2) # mz assuming charge of 2

    return([flat_peptide_abundances, flat_organism_profiles,
            organism_index, peptide_mass, list(retention_times),
            flat_organism_protein_pools])

def check_pool_or_tax_ambiguous(peptide_index_vals, gen_dataset):
    '''
    Some peptides are observed more than once because they are from the same organism,
    but in different protein pools. Some are observed twice because they are from different organisms.
    This function goes into the peptide indexes which are duplicates (or more than duplicates) and
    then alters the description of the organism index or the protein pool index to be 'ambigious'.
    '''
    all_taxa = []
    all_protein_pool = []
    # subset the dataset based on the peptide indexes
    for i in peptide_index_vals:
    # check the number of unique protein pools
    # check the number of unique taxa
        all_taxa.append(gen_dataset[2][i])
        all_protein_pool.append(gen_dataset[5][i])

    length_unique_taxa = len(set(all_taxa))
    length_unique_pools = len(set(all_protein_pool))

    if length_unique_taxa > 1 and length_unique_pools > 1:
        return('both ambiguous')
    if length_unique_taxa > 1 and length_unique_pools == 1:
        return('org ambiguous')
    if length_unique_taxa == 1 and length_unique_pools > 1:
        return('pool ambiguous')

def get_observed_dataset(gen_dataset_in):
    '''
    There will be duplicate peptides observed, this function makes an 'observed dataset' where
    those peptides are duplicated
    '''
    # make a copy (this is from Noor Youssef, Aug. 2019)
    gen_dataset = [x[:] for x in gen_dataset_in]

    for peptide in gen_dataset[1]:
        # is this peptide only found once?
        number_instances = gen_dataset[1].count(peptide)
#         print(peptide)

        if number_instances == 1:
            continue

        # get the indexes of all peptides
        peptide_index = [i for i, e in enumerate(gen_dataset[1]) if e == peptide]

        # get the peptide intensity of the first value (0)
        new_intensity_value = gen_dataset[0][peptide_index[0]]

        ## input function here above(check_pool_or_tax_ambiguous)
        which_to_be_amb = check_pool_or_tax_ambiguous(peptide_index_vals = peptide_index,
                                                      gen_dataset = gen_dataset)
        # designate which component is ambiguous
        if which_to_be_amb == 'org ambiguous':
            gen_dataset[2][peptide_index[0]] = 'ambiguous'
        if which_to_be_amb == 'pool ambiguous':
            gen_dataset[5][peptide_index[0]] = 'ambiguous'
        if which_to_be_amb == 'both ambiguous':
            gen_dataset[2][peptide_index[0]] = 'ambiguous'
            gen_dataset[5][peptide_index[0]] = 'ambiguous'

        # start the index to be removed from 1, not 0, so it keeps the original
        for to_be_removed_index in reversed(peptide_index):
            # if it's at the first element, then don't add it again
            if to_be_removed_index == peptide_index[0]:
                continue
            # sum the intensity values for the peptides that are the same
            new_intensity_value += gen_dataset[0][to_be_removed_index]

            # remove this redundant peptide from every other element of the list of lists
            for i in range(len(gen_dataset)):
                del gen_dataset[i][to_be_removed_index]

        gen_dataset[0][peptide_index[0]] = new_intensity_value

    return(gen_dataset)

def write_pd_dataset(n_orgs, n_seqs,
                  number_of_seqs_vec = [250, 1000, 2500, 10000],
                  rt_low = 0,
                  rt_upp = 120):

    dataset_true = write_dataset(n_orgs = n_orgs, n_seqs = n_seqs,
                  number_of_seqs_vec = number_of_seqs_vec,
                  rt_low = rt_low,
                  rt_upp = rt_upp)

    dataset_obs = get_observed_dataset(dataset_true)

    pd_dataset_true = DataFrame(dataset_true)
    pd_dataset_obs = DataFrame(dataset_obs)

    pd_dataset_trans_obs = pd_dataset_obs.transpose()
    pd_dataset_trans_true = pd_dataset_true.transpose()

    pd_dataset_trans_obs.columns = ['ion_intensity', 'peptide_sequence',
                                    'organism', 'mz', 'rts', 'protein_pool']
    pd_dataset_trans_true.columns = ['ion_intensity', 'peptide_sequence',
                                    'organism', 'mz', 'rts', 'protein_pool']
    return([pd_dataset_trans_obs,
            pd_dataset_trans_true])

##### GENERATIVE MODEL (part 2)

### output should be another pandas dataframe that is of the peptides identified and their corresponding
### masses

def _expire_ion(ms_list, time_span, ms_run_time):
    ion_injection_time = ms_list[2]
    ion_bool = (ion_injection_time + time_span) > ms_run_time
    return ion_bool

# filtering function for selecting ion parcels
def _filter_ion_parcel(peptide_df, filtered_injection_bin_i):
    if type(filtered_injection_bin_i) != pd._libs.interval.Interval:
        raise NameError('Error with retention time binning. Check pandas version for interval class support (>0.2.0).')
    ion_parcel = peptide_df[(peptide_df.rt_upper > filtered_injection_bin_i.left) & 
                            (peptide_df.rt_upper <= filtered_injection_bin_i.right) |
        (peptide_df.rt_lower > filtered_injection_bin_i.left) &
        (peptide_df.rt_lower <= filtered_injection_bin_i.right) |
        (peptide_df.rt_lower <= filtered_injection_bin_i.left) &
        (peptide_df.rt_upper >= filtered_injection_bin_i.right)]
    return(ion_parcel)

def format_injection_bins(peptide_unique,
                         max_injection_time = 5, ion_peak_width = 1):
    '''
    the peptide output above needs to be formatted with 'injection bins', replicating the injection
    mechanism of mass spectra
    '''
    injection_bins_ranges = np.arange(start = peptide_unique['rts'].min(),
                           stop = peptide_unique['rts'].max() + peptide_unique['rts'].max()*0.05,
                           step = max_injection_time)
    injection_bins = pd.cut(peptide_unique['rts'].tolist(),
        injection_bins_ranges,
        right = True,
        include_lowest = True)
    peptide_unique.loc[:,'injection_bins'] = pd.Series(injection_bins,
                                                       index = peptide_unique.index)
    #%% Range for ion peak width.

    peptide_unique = peptide_unique.assign(rt_upper = pd.Series(peptide_unique['rts'] + ion_peak_width), 
                                           rt_lower = pd.Series(peptide_unique['rts'] - ion_peak_width))

    return(peptide_unique)

def run_stoch_sim(peptide_unique, ion_peak_width = 1,
                  top_n_ions = 12, 
                  precursor_selection_window = 3, time_span = 30):
    '''
    mass spec sampling model. this model simulates data dependent acquisition of mass spectra using
    a top n method. it includes dynamic exclusion.
    '''

    peptide_unique_sorted = peptide_unique.sort_values(['rts'])

    sorted_unique_injection_bins = peptide_unique_sorted['injection_bins'].unique()

    # Dynamic exclusion list of ions. Numbers are upperbound, lowerbound, and time.
    dynamic_exclusion_ions = [[0, 0, 0]]

    # Counter variable, which increases along the for loop.
    ms_run_time = 0 # minutes.

    # Initialize this cofragmentation counter column, which is added to iteratively.
    # Initialize the simulation tracker variables
    total_dynamic = []
    number_ions_selected = []


    # in an MS experiment of top N, a precursor ion scan happens once and then a
    # product ion scan happens N times.
    # This line subsets the injection bins
    subset_sorted_unique_injection_bins = sorted_unique_injection_bins[::top_n_ions]

    sampled_ion_parcels = []

    #%% Looping through ion parcels (via injection_bins). Subsamples based on ions in parcel.
    for injection_bin_i in range(0, len(set(subset_sorted_unique_injection_bins))):

        filtered_injection_bin = subset_sorted_unique_injection_bins[injection_bin_i]

        # Ion parcel includes any ion which has an overlapping retention time bin (referred to with
        # rt_upper/rt_lower columns) within the injection_bin.
        ion_parcel = _filter_ion_parcel(peptide_df = peptide_unique,
                                           filtered_injection_bin_i = filtered_injection_bin)

        # Filter out rows within a certain bins from dynamic exclusion.
        # Goes through each ion in ion_parcel and removes peptides that are dynamically excluded.
        for row in dynamic_exclusion_ions:
            lower_bound = row[0]
            upper_bound = row[1]
            ion_parcel = ion_parcel[(ion_parcel.mz < lower_bound) | (ion_parcel.mz > upper_bound)]

        # Sometimes all ions will be blocked out by the dynamic exclusion, if so, continue.
        if(ion_parcel.shape[0] == 0):
            continue

        # Sample ions
            # If there are fewer peptides than top_n, take # peptides.
        if ion_parcel.shape[0] < top_n_ions:
            top_n_sample = ion_parcel.shape[0]
            #sampled_ion_parcel = ion_parcel.sample(top_n_sample)
            sampled_ion_parcel = ion_parcel
        else:
            # select the top most intense ions
            sampled_ion_parcel = ion_parcel.sort_values(by = 'ion_intensity',
                                                        ascending = False)[:top_n_ions]

        # appending the sampled ions to a list
        sampled_ion_parcels.append(sampled_ion_parcel)

        # Updating ms run time
        ms_run_time = filtered_injection_bin.right

        # Updating the dynamic exclusion list based on sampled_ions.
        for row_i in range(0, sampled_ion_parcel.shape[0]):
            dynamic_exclusion_row = [sampled_ion_parcel['mz'].values[row_i] - precursor_selection_window/2, sampled_ion_parcel['mz'].values[row_i] + precursor_selection_window/2, sampled_ion_parcel['rts'].values[row_i]]
            dynamic_exclusion_ions.append(dynamic_exclusion_row)

        # Removing ions which have expired after some time_span.
        dynamic_exclusion_ions[:] = [tup for tup in dynamic_exclusion_ions if _expire_ion(tup,
                                                                                          time_span = time_span,
                                                                                          ms_run_time = ms_run_time)]

    # format all sampled ions into a pd df
    formatted_ions = pd.concat(sampled_ion_parcels)

    return(formatted_ions)

def run_gen_model(n_orgs = 10, n_seqs = 200,
                  number_of_seqs_vec = [250, 1000, 2500, 10000],
                  rt_low = 0,
                  rt_upp = 120, max_injection_time = 0.00833333, ion_peak_width = 0.5,
                 top_n_ions = 12, precursor_selection_window = 3, time_span = 30):

    print('Running generative model...')
    gen_dataset = write_pd_dataset(n_orgs = n_orgs,
                     n_seqs = n_seqs,
                     number_of_seqs_vec = number_of_seqs_vec,
                     rt_low = rt_low,
                     rt_upp = rt_upp)
    format_gen_dataset = format_injection_bins(peptide_unique=gen_dataset[0],
                         max_injection_time = max_injection_time,
                         ion_peak_width = ion_peak_width)

    print('Running sampling model...')
    ms_subsample = run_stoch_sim(peptide_unique = format_gen_dataset,
                                  ion_peak_width = ion_peak_width,
                                  top_n_ions = top_n_ions,
                                  precursor_selection_window = precursor_selection_window,
                                  time_span = time_span)

    return([ms_subsample, gen_dataset[0], gen_dataset[1]])

### Estimating true proteome from simulated metaproteome:

def get_multiple_rmse(n_trials, file_name, n_orgs = 100, n_seqs = 2000,
                  number_of_seqs_vec = [100000, 200000, 500000, 1000000],
                  rt_low = 0, rt_upp = 90, max_injection_time = 0.00833333,
                  ion_peak_width = 0.5, top_n_ions = 12,
                  precursor_selection_window = 3, time_span = 0.5):

    sampled_out = []
    observed_true = []
    actual_true = []

    for trial in range(n_trials):
        print('Trial: ', trial, flush = True)
        gen_model = run_gen_model(n_orgs, n_seqs,
                      number_of_seqs_vec,
                      rt_low, rt_upp, max_injection_time,
                      ion_peak_width, top_n_ions,
                      precursor_selection_window, time_span)

        gen_model[0]['trial'] = str(trial)
        gen_model[1]['trial'] = str(trial)
        gen_model[2]['trial'] = str(trial)

        sampled_out.append(gen_model[0])
        observed_true.append(gen_model[1])
        actual_true.append(gen_model[2])


    sampled_models_out = pd.concat(sampled_out)
    observed_true_out = pd.concat(observed_true)
    actual_true_out = pd.concat(actual_true)

    sampled_models_out.to_csv('sampled_' + file_name)
    observed_true_out.to_csv('observed_' + file_name)
    actual_true_out.to_csv('actual_' + file_name)

get_multiple_rmse(n_trials = 15, file_name = 'proper_observed_orgs30_n_seqs2000_5prot_groups_updated_ambiguous.csv', n_orgs = 30, n_seqs = 2000, number_of_seqs_vec = [15000, 50000, 100000, 250000, 500000], rt_upp = 260)

