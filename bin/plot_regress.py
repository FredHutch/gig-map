#!/usr/bin/env python3

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import plotly.express as px


def main():

    df = pd.read_csv("regress.results.csv")

    # Calculate the fold change relative to the intercept
    df = (
        pd.concat([
            calc_fold_change(feature_df)
            for _, feature_df in df.groupby("feature")
        ])
        .query("parameter != '(Intercept)'")
        .pipe(add_columns)
    )

    # Make some plots
    for param, param_df in df.groupby("parameter"):
        plot_scatter(param_df, "fold_change", "neg_log10_pvalue", f"{param}.fold_change_vs_pvalue.scatter")
        plot_scatter(param_df, "fold_change", "neg_log10_qvalue", f"{param}.fold_change_vs_qvalue.scatter")

    # Save the data that was used for plotting
    df.sort_values(by="pvalue").to_csv("association.csv", index=None)


def plot_scatter(df: pd.DataFrame, x: str, y: str, output_prefix: str):
    fig = px.scatter(
        data_frame=df,
        x=x,
        y=y,
        template="simple_white",
        hover_name="feature",
        hover_data=["pvalue", "qvalue", "Estimate", "fold_change"],
        labels=dict(
            Estimate="Estimated Coefficient of Association",
            fold_change="Fold Change",
            pvalue="p-value",
            qvalue="q-value",
            neg_log10_pvalue="-log10(pvalue)",
            neg_log10_qvalue="-log10(qvalue)",
        )
    )
    fig.write_html(output_prefix + ".html")
    fig.write_image(output_prefix + ".png")


def add_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.assign(
        qvalue=lambda d: multipletests(d["pvalue"], 0.01, "fdr_bh")[1],
        neg_log10_pvalue=lambda d: d["pvalue"].apply(lambda v: -np.log10(v)),
        neg_log10_qvalue=lambda d: d["qvalue"].apply(lambda v: -np.log10(v)),
    )


def calc_fold_change(bin_df: pd.DataFrame) -> pd.DataFrame:
    # Get the estimate for each parameter
    estimate = bin_df.set_index("parameter")["Estimate"]

    # The fold change is the estimate as a proportion of the maximum value
    intercept = estimate.loc["(Intercept)"]
    fold_change = estimate.apply(
        lambda v: v / np.max([intercept, intercept + v])
    )
    return bin_df.assign(fold_change=bin_df["parameter"].apply(fold_change.get))


if __name__ == "__main__":
    main()
