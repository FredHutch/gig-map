#!/usr/bin/env python3

# Import libraries
import gzip
import logging
import pandas as pd
from sklearn.metrics import pairwise_distances
import sys
from time import time
from uuid import uuid4


class Timer:
    """Log the amount of time taken to perform a task."""

    def __init__(self, logger=None):

        # Record the time that the timer was created
        self.start_time = time()

        # Attach the logger, if any
        self.logger = logger

    def stop(self, msg):
        """Report the amount of time elapsed."""

        # Format the message with the amount of time elapsed
        msg = f"{msg} - {self.elapsed():.2E} seconds elapsed"

        # If there is a logger
        if self.logger is not None:

            # Log the message
            self.logger.info(msg)

        # If there is no logger
        else:

            # Print the message
            print(msg)

    def elapsed(self):
        """Return the number of seconds since the timer started."""
        return time() - self.start_time


def calc_dm(df, metric='cosine'):
    """Return a square DataFrame of pairwise distances between rows."""

    return pd.DataFrame(
        pairwise_distances(
            df.values,
            metric=metric,
            # Use all available cores
            n_jobs=-1
        ),
        index=df.index.values,
        columns=df.index.values
    )


def order_rows(
    # Table with values used to sort the rows
    df:pd.DataFrame,
    # Distance metric
    metric="cosine",
    # Logging instance
    logger=None,
    # Optionally allow the user to specify a row to start the order
    starting_row=None,
    # Number of seconds between logging messages
    log_interval=10
):
    """
    Order the rows of a DataFrame by greedily selecting adjacencies.
    """
    
    # Compute the distance matrix
    t = Timer(logger=logger)
    dm = calc_dm(df, metric=metric)
    t.stop(f"Calculated pairwise {metric} distances for {dm.shape[0]} rows")

    # Sort the distances for each row
    t = Timer(logger=logger)
    dists = {
        row_name: row_dists.sort_values()
        for row_name, row_dists in dm.iterrows()
    }
    t.stop(f"Sorted distances for all rows")

    # If the user did not specify a starting row
    if starting_row is None:

        # Set the starting row as the one with the largest sum
        starting_row = df.mean(axis=1).sort_values().index.values[-1]
        logger.info(f"Order will start with row {starting_row}")

    # If the user did specify a starting row
    else:

        # Make sure that the row is in the table
        assert starting_row in df.index.values, f"Row {starting_row} not found"

    # Set up the row order
    row_order = [starting_row]
    seen = set(row_order)

    t = Timer(logger=logger)
    
    while len(seen) < len(dists):

        # Find the next row to add to the ordered list
        next_row = best_match(
            row_order[-1],
            dists,
            seen
        )

        # Add that row to the list
        row_order.append(next_row)
        seen.add(next_row)
            
        if t.elapsed() > log_interval:
            t.stop(f"Ordered {len(seen):,} genes")
            t = Timer(logger=logger)

    # Return the final list
    return row_order


def best_match(query_row, dists, seen):
    """Return the most similar row which has not yet been seen."""
    
    for row_name, dist in dists[query_row].items():
        
        if row_name not in seen:
            
            return row_name


# If this is being run as a script
if __name__ == "__main__":

    ##################
    # SET UP LOGGING #
    ##################

    # Set the level of the logger to INFO
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [Order Genes] %(message)s'
    )
    logger = logging.getLogger(str(uuid4()))
    logger.setLevel(logging.INFO)

    # Write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    # Get the user-supplied positional arguments
    input_fp = sys.argv[1]
    logger.info(f"Reading input from {input_fp}")
    output_fp = sys.argv[2]

    # Read in all of the alignment information
    df = pd.read_csv(
        input_fp
    # Pivot to wide format
    ).pivot_table(
        columns="genome",
        index="sseqid",
        values="pident"
    ).fillna(
        0
    )

    # Get the order of rows based on linkage clustering
    t = Timer(logger=logger)
    gene_order = order_rows(
        df,
        logger=logger
    )
    t.stop("Overall gene ordering")

    logger.info(f"Writing output to {output_fp}")

    # Write out to a file
    with gzip.open(
        output_fp,
        "wt"
    ) as handle:
        handle.write("\n".join(gene_order))

    logger.info("DONE")