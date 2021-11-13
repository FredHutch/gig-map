#!/usr/bin/env python3

import shutil
import urllib.request as request
from urllib.error import URLError
from contextlib import closing


remote_path = "${ftp_url}"
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
    if "${params.skip_missing_ftp}" == "true":
        print("Missing FTP paths will be ignored")
    else:
        print("Raising error")
        raise e
