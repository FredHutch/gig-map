#!/usr/bin/env python3

import gzip


class Shard:

    def __init__(self):

        self.shard_size = ${params.map_batchsize}
        print(f"Breaking up queries into batches of {self.shard_size:,} genes each")
        self.shard_ix = 0
        self.handle = None
        self.open_handle()
        self.gene_ix = 0

    def open_handle(self):
        self.handle = gzip.open(f"queries.shard.{self.shard_ix}.fasta.gz", "wt")

    def increment(self):
        self.handle.close()
        self.shard_ix += 1
        self.open_handle()

    def write(self, line):
        if line[0] == ">":
            self.gene_ix += 1
            if self.gene_ix % self.shard_size == 0:
                self.increment()
        self.handle.write(line)

    def read(self):
        with gzip.open("queries.fasta.gz", "rt") as handle:
            for line in handle:
                self.write(line)


shard = Shard()
shard.read()
