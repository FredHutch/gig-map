#!/usr/bin/env python3

# Import libraries
import gzip
import logging
import pandas as pd
from random import choice
from scipy.spatial import distance
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
        msg = f"{msg} - {time() - self.start_time:.2E} seconds elapsed"

        # If there is a logger
        if self.logger is not None:

            # Log the message
            self.logger.info(msg)

        # If there is no logger
        else:

            # Print the message
            print(msg)


def calc_dm(df, metric='dice'):
    """Return a square DataFrame of pairwise distances between rows."""

    return pd.DataFrame(
        distance.squareform(
            distance.pdist(
                df.values,
                metric=metric
            )
        ),
        index=df.index.values,
        columns=df.index.values
    )


def order_rows(
    df:pd.DataFrame,
    metric="dice",
    logger=None,
    n_iter=5,
):
    """
    Order the rows of a DataFrame by greedily selecting adjacencies.
    """
    
    # Compute the distance matrix
    t = Timer(logger=logger)
    dm = calc_dm(df, metric=metric)
    t.stop(f"Calculated pairwise distances for {dm.shape[0]} genes")

    # Keep track of the best score and the best order
    best_score, best_order = None, None

    # Iteratively sort the distance matrix
    for iter_i in range(n_iter):
        
        # In each iteration, sort the rows by greedy_linearization
        t = Timer(logger=logger)
        new_score, new_order = greedy_linearization(dm)
        t.stop(f"Iteration {iter_i}: Score = {new_score}")
        
        # If the score is better
        if best_score is None or new_score < best_score:
            
            # Save it as the best
            best_score, best_order = new_score, new_order

    logger.info(f"Top score: {best_score}")

    # After all iterations are done, return the best order that was found
    return best_order


def greedy_linearization(dm):
    """Return the order of rows by greedily selecting adjacent rows."""
    
    # Make a list, starting with a random row
    row_order = [
        choice(dm.index.values)
    ]
    
    # Keep track of the total score
    total_score = 0
    
    # While there are more rows to add to the list
    while len(row_order) < dm.shape[0]:

        # Get the candidates of rows to choose from
        candidates = dm.loc[
            row_order[-1]
        ].drop(
            index=row_order
        )
        
        # Get the best match for the final row
        best_match = candidates.idxmin()
        
        # Add to the score
        total_score += candidates.min()
        
        # Add to the order
        row_order.append(best_match)
    
    return total_score, row_order


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
    # Convert to a presence / absence boolean
    ).applymap(
        lambda v: int(v > 0)
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