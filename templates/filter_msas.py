#!/usr/bin/env python3

import gzip
import pandas as pd
from pathlib import Path
import logging

logger = logging.getLogger('gig-map')
logger.setLevel(logging.INFO)
# Write to STDOUT
consoleHandler = logging.StreamHandler()
logger.addHandler(consoleHandler)


def read_gzipped_fasta(fasta_path: str) -> dict:
    sequences = {}
    with gzip.open(fasta_path, 'rt') as f:
        seq_id = None
        seq_chunks = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seq_id is not None:
                    sequences[seq_id] = ''.join(seq_chunks)
                seq_id = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if seq_id is not None:
            sequences[seq_id] = ''.join(seq_chunks)
    return sequences


def filter_msa(msa_fp: Path, min_prop_sites: float, min_prop_genomes: float, max_raxml_sites: int) -> pd.DataFrame:
    """
    Filter the MSA to only include sites which meet the minimum proportion of sites,
    and only include genomes which have at least this proportion of the sites present.
    """
    logger.info(f"Filtering MSA {msa_fp} with min_prop_sites={min_prop_sites}, min_prop_genomes={min_prop_genomes}")

    # Read in the MSA from gzipped FASTA
    seqs = read_gzipped_fasta(msa_fp)
    msa = pd.DataFrame([list(seq) for seq in seqs.values()], index=list(seqs.keys()))
    logger.info(f"Read in MSA with {msa.shape[0]} genomes and {msa.shape[1]} sites from gzipped FASTA")
    
    # Iteratively filter sites and genomes until convergence
    prev_shape = (0, 0)
    msa_filtered = msa.copy()
    while msa_filtered.shape != prev_shape:
        prev_shape = msa_filtered.shape
        # Filter genomes based on min_prop_genomes
        genome_non_missing = (msa_filtered != 'N').sum(axis=1)
        genome_prop_non_missing = genome_non_missing / msa_filtered.shape[1]
        filtered_genomes = genome_prop_non_missing[genome_prop_non_missing >= min_prop_genomes].index
        msa_filtered = msa_filtered.loc[filtered_genomes]
        logger.info(f"Filtered MSA to {msa_filtered.shape[0]} genomes after applying min_prop_genomes")
        # Filter sites based on min_prop_sites
        site_non_missing = (msa_filtered != 'N').sum(axis=0)
        site_prop_non_missing = site_non_missing / msa_filtered.shape[0]
        filtered_sites = site_prop_non_missing[site_prop_non_missing >= min_prop_sites].index
        msa_filtered = msa_filtered[filtered_sites]
        logger.info(f"Filtered MSA to {msa_filtered.shape[1]} sites after applying min_prop_sites")

    # If there are more than max_raxml_sites sites
    if msa_filtered.shape[1] > max_raxml_sites:
        logger.info(f"Truncating down to {max_raxml_sites:,} sites")
        msa_filtered = msa_filtered.iloc[:, :max_raxml_sites]

    return msa_filtered


def main():

    # Parameters from config
    min_prop_sites = ${params.raxml_min_prop_sites}
    min_prop_genomes = ${params.raxml_min_prop_genomes}
    min_raxml_genomes = ${params.min_raxml_genomes}
    max_raxml_sites = ${params.max_raxml_sites}
    
    # Input MSA file
    input_msa_fp = Path("input").rglob("*.msa.gz").__next__()
    
    # Output filtered MSA file
    output_msa_fp = input_msa_fp.name
    
    # Filter the MSA
    filtered_msa = filter_msa(input_msa_fp, min_prop_sites, min_prop_genomes, max_raxml_sites)

    if filtered_msa.shape[0] == 0 or filtered_msa.shape[1] == 0:
        logger.info("No genomes or sites remain after filtering. Exiting.")
        return
    
    if filtered_msa.shape[0] < min_raxml_genomes:
        logger.info(f"The number of genomes passing the filter ({filtered_msa.shape[0]}) does not meet the minimum threshold ({min_raxml_genomes})")
        return
    
    # Save the filtered MSA in FASTA format
    with gzip.open(output_msa_fp, 'wt') as fasta_out:
        for genome, row in filtered_msa.iterrows():
            sequence = ''.join(row.tolist())
            fasta_out.write(f'>{genome}\\n{sequence}\\n')

    logger.info(f"Saved filtered MSA in FASTA format to {output_msa_fp}")


if __name__ == "__main__":
    main()