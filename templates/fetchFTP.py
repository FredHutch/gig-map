#!/usr/bin/env python3

import gzip
import shutil
import urllib.request as request
from urllib.error import URLError
from contextlib import closing

# The path to download will be filled in by Nextflow
remote_path = "${ftp_url}"

# We will make a backup path which references the RefSeq record instead of GenBank
backup_path = remote_path.replace("GCA", "GCF")


def download_path(remote_path:str, skip_missing_ftp=True):
    local_path = remote_path.rsplit("/", 1)[-1]
    print(f"Downloading from {remote_path}")
    print(f"Local path is {local_path}")

    try:
        with closing(
            request.urlopen(
                remote_path,
                timeout=60
            )
        ) as r:
            with open(local_path, 'wb') as f:
                shutil.copyfileobj(r, f)
    except URLError as e:
        print("The URL appears to not be valid")
        if skip_missing_ftp:
            print("Missing FTP paths will be ignored")
            return False
        else:
            print("Raising error")
            raise e

    # If the file is expected to be gzip compressed
    if local_path.endswith('.gz'):

        # Try to open the file
        with gzip.GzipFile(local_path, mode="r") as handle:

            # Try to read it
            handle.read()

    return True


# If we can't download the first option
if not download_path(remote_path, skip_missing_ftp=False):

    # Try the second option
    print(f"Trying backup path: {backup_path}")
    download_path(backup_path, skip_missing_ftp="${params.skip_missing_ftp}" == "true")
