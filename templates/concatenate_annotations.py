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
    # and yield a series of:
    # ['key1', 'key2'], 1
    for path, val in flatten_dict(dat):
        # First see if ['key2'] can be used for ['key1', 'key2']
        # If not, keep expanding until we get a match
        kw = None
        for i in range(1, len(path)):
            kw = "_".join(path[-i:])
            if key_map.get(kw) == path:
                break
            elif kw not in key_map:
                key_map[kw] = path
                break
        if isinstance(val, str) and "\\n" in val:
            continue
        output[kw] = val

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
                    for x, y in flatten_dict(v, path=path + [str(i), kw]):
                        yield x, y
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
