#!/usr/bin/env python3
import pandas as pd
import json
import os

# Set the name of the output file based on the output_prefix value (genes/genomes)
output_fp = "${output_prefix}.annot.csv.gz"

key_map = dict()


def parse_jsonl(fp):
    with open(fp) as handle:
        dat = json.load(handle)
    # The data returned is in a nested dict that we want to turn into a wide DataFrame
    output = {}

    # flatten_dict will take:
    # {'key1': {'key2': 1}}
    # and return an iterator of (e.g.):
    # ['key1', 'key2'], 1
    for path, val in flatten_dict(dat):
        # Add the value to the output
        output["_".join(path)] = val

    # Optionally merge the accession and assemblyName fields
    assert "accession" in output, list(output.keys())
    assert "assemblyInfo_assemblyName" in output, list(output.keys())
    output["genome_id"] = "_".join([
        output["accession"],
        output["assemblyInfo_assemblyName"],
        "genomic.fna.gz"
    ])

    return output


def flatten_dict(dat: dict, path=None):
    if path is None:
        path = []
    for kw, val in dat.items():
        if isinstance(val, dict):
            for x, y in flatten_dict(val, path=path + [kw]):
                yield x, y
        elif isinstance(val, list):
            for i, v in enumerate(val):
                if isinstance(v, dict):
                    for x, y in flatten_dict(v, path=path + [kw, str(i)]):
                        yield x, y
                else:
                    yield path + [str(i)], v
        else:
            yield path + [kw], val


# If there are any annotations available
if os.path.exists('annotations'):

    # Read in all of the inputs
    df = pd.concat((
        [
            pd.read_csv(
                os.path.join('annotations', fp)
            )
            for fp in os.listdir('annotations')
            if fp.endswith('.csv.gz')
        ]
        + [pd.DataFrame([
            parse_jsonl(
                os.path.join('annotations', fp)
            )
            for fp in os.listdir('annotations')
            if fp.endswith('.jsonl')
        ])]
    ))

    # Write out the table
    df.to_csv(
        output_fp,
        index=None
    )

# If there are no annotations available
else:

    # Write out an empty file
    with open(output_fp, 'wt') as handle:
        handle.write("${output_prefix}_id\\n")
