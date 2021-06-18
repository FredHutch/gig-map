#!/usr/bin/env python3
"""Script to reformat gene annotations from geneshot.results.hdf5 to CSV as needed for gig-map."""

import numpy as np
import pandas as pd
import argparse

###################
# PARSE ARGUMENTS #
###################

# Create the parser
parser = argparse.ArgumentParser(
    description="Reformat gene annotations from geneshot.results.hdf5 to CSV as needed for gig-map"
)

# Add the arguments
parser.add_argument(
    '--input',
    type=str,
    required=True,
    help='Geneshot results in HDF5 format (.results.hdf5)'
)
parser.add_argument(
    '--output',
    type=str,
    required=True,
    help='Path to write annotations in CSV format (.csv.gz)'
)

# Parse the arguments
args = parser.parse_args()

# Make sure that the output path has the expected extension
assert args.output.endswith(".csv.gz"), "--output must end with .csv.gz"


def format_gene_annot(r):
    """Function to format a single row of the gene annotation table."""

    # Set up a dict with the output
    output = dict()

    # Set the gene ID
    output['gene_id'] = r['gene']

    # Annotate a couple fields directly
    for k in ['CAG', 'length']:
        output[k] = r[k]

    # If there is a taxonomic annotation
    if 'tax_name' in r.index.values:

        # Which is not null
        if not pd.isnull(r['tax_name']):

            # Add the name by itself
            output['taxonomic_assignment'] = r['tax_name']

    # If there is a functional annotation
    if 'eggNOG_desc' in r.index.values:

        # Which is not null
        if not pd.isnull(r['eggNOG_desc']):

            # Add the name by itself
            output['functional_assignment'] = r['eggNOG_desc']

    # Finally, format a combined name
    output['combined_name'] = format_combined_name(
        output,
        check_for_function='eggNOG_desc' in r.index.values
    )

    # Return a pandas.Series based on the output
    return pd.Series(output)


def format_combined_name(output, check_for_function=False):
    """Format a single human-readable name for a gene."""

    # If the dataset overall contains functional annotations
    if check_for_function:

        # The combined name will start with the functional annotation
        combined_name = output.get(
            'functional_assignment',
            # Defaulting to the string below for genes which lack those annotations
            "Unannotated sequence"
        )

        # Now append the unique gene ID
        combined_name = f"{combined_name} (ID: {output['gene_id']})"

    # Otherwise, if the dataset did not have any functional annotation performed
    else:

        # The combined name will start with the gene ID
        combined_name = output["gene_id"]

    # If there is a taxonomic annotation
    if output.get('taxonomic_assignment') is not None:

        # Add the name of the organism
        combined_name = f"{combined_name} [{output['taxonomic_assignment']}]"

    # End with the gene length
    combined_name = f"{combined_name} {output['length']}aa"

    # Return the name
    return combined_name


# Open a connection to the HDF store
with pd.HDFStore(args.input, "r") as store:

    # Read in the table with all gene annotations
    gene_annots = pd.read_hdf(
        store,
        "/annot/gene/all"
    )

    # If there are eggNOG annotations
    if "eggNOG_desc_ix" in gene_annots.columns.values:

        # Read in the mapping of eggNOG_desc_ix to human-readable names
        eggnog_name_dict = pd.read_hdf(
            store,
            "/ref/eggnog_desc"
        ).set_index(
            "eggnog_desc_ix"
        )["eggnog_desc"].to_dict()

        # Add the human-readable string to the DataFrame
        gene_annots = gene_annots.assign(
            eggNOG_desc = gene_annots.eggNOG_desc_ix.apply(
                lambda v: eggnog_name_dict.get(int(v)) if not pd.isnull(v) else None
            )
        )

    # Format the output based on the content of per-gene annotations
    gene_annots = gene_annots.apply(
        format_gene_annot,
        axis=1
    ).set_index(
        "gene_id"
    ).sample(3)

    # If there are results of association tests at the CAG-level with experimental design
    if "/stats/cag/corncob" in store:

        # Read in the stats table
        stats_df = pd.read_hdf(
            store,
            "/stats/cag/corncob"
        ).query(
            "parameter != '(Intercept)'"
        )

        # Clean up the columns which we'll output
        stats_df = stats_df.assign(
            neg_log10_qvalue = -1 * stats_df['q_value'].apply(np.log10)
        ).drop(
            columns=['std_error', 'p_value', 'q_value']
        )

        # Add the paramter name to the columns
        stats_df = pd.concat([
            param_df.rename(
                columns=dict(cag='CAG')
            ).set_index(
                'CAG'
            ).drop(
                columns=['parameter']
            ).rename(
                columns=lambda col_name: f"{param} - {col_name}"
            )
            for param, param_df in stats_df.groupby("parameter")
        ], axis=1)

        # Add the CAG stats results to the gene annotation table
        gene_annots = pd.merge(gene_annots, stats_df, on='CAG')

    # Write out to a file
    gene_annots.to_csv(args.output)

    print(f"Wrote out to {args.output}")
