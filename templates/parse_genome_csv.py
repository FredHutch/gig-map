#!/usr/bin/env python3

import pandas as pd
import re

df = pd.read_csv("input.csv")

# Make a list of URLs to download
url_list = []

# Function to format the path to the genome FASTA
def format_ftp(ftp_prefix):

    assert isinstance(ftp_prefix, str)
    assert ftp_prefix.startswith("ftp://")
    assert "/" in ftp_prefix
    assert not ftp_prefix.endswith("/")

    # The ID of the assembly is the final directory name
    id_str = ftp_prefix.rsplit("/", 1)[-1]

    # Return the path to the genome FASTA
    return f"{ftp_prefix}/{id_str}${params.parse_genome_csv_suffix}"

# Populate a list with formatted annotations for each file
# which will be downloaded from NCBI
annotation_list = []

# Function to format the annotation for each row
def format_annotation(r):

    # Get the name of the file which will be downloaded
    annots = dict(
        genome_id=r["GenBank FTP"].rsplit("/", 1)[-1] + "${params.parse_genome_csv_suffix}"
    )

    # Rename fields
    for k, n in [
        ("#Organism Name", "Organism"),
        ("Strain", "Strain"),
        ("BioSample", "BioSample"),
        ("BioProject", "BioProject"),
        ("Assembly", "Assembly"),
        ("Size(Mb)", "Size_Mb"),
        ("GC%", "GC_percent"),
        ("CDS", "CDS")
    ]:
        annots[n] = r.get(k)

    # Format a combined name
    combined_name = annots["Organism"]
    if isinstance(annots["Strain"], str):
        combined_name = combined_name + "(" + annots["Strain"] + ")"
    annots["Formatted Name"] = combined_name + "[" + annots["Assembly"] + "]"

    return annots

# Iterate over each row in the table
for _, r in df.iterrows():
    
    # If there is no value in the 'GenBank FTP' column
    if pd.isnull(r['GenBank FTP']):
    
        # Skip it
        continue
    
    # If the 'GenBank FTP' column doesn't start with 'ftp://'
    elif not r['GenBank FTP'].startswith('ftp://'):
    
        # Skip it
        continue

    # Otherwise
    else:

        # Format the path to the genome in that folder
        url_list.append(
            format_ftp(r['GenBank FTP'])
        )

        # Format the annotations
        annotation_list.append(
            format_annotation(r)
        )

# Write the list to a file
with open("url_list.txt", "w") as handle:

    # Each URL on its own line
    handle.write("\\n".join(url_list))

# Write the annotations to a CSV
pd.DataFrame(annotation_list).set_index("genome_id").to_csv(
    "genome_annotations.csv.gz"
)