#!/usr/bin/env python3

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import plotly.express as px


def main():

    df = (
        pd.read_csv("radEmu.results.csv")
        .rename(
            columns=dict(
                estimate="Estimate",
                pval="pvalue",
                covariate="parameter",
                category="feature"
            )
        )
        .pipe(add_columns)
    )

    # Make some plots
    for param, param_df in df.groupby("parameter"):
        plot_scatter(param_df, "Estimate", "neg_log10_pvalue", f"{param}.estimate_vs_pvalue.scatter")
        plot_scatter(param_df, "Estimate", "neg_log10_qvalue", f"{param}.estimate_vs_qvalue.scatter")

    # Save the data that was used for plotting
    df.sort_values(by="pvalue").to_csv("association.csv", index=None)


def plot_scatter(df: pd.DataFrame, x: str, y: str, output_prefix: str):
    fig = px.scatter(
        data_frame=df,
        x=x,
        y=y,
        template="simple_white",
        hover_name="feature",
        hover_data=["pvalue", "qvalue", "Estimate"],
        labels=dict(
            Estimate="Estimated Coefficient of Association",
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


if __name__ == "__main__":
    main()
