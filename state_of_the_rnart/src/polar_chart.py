import os
from typing import List

import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.preprocessing import MinMaxScaler

ALL_MODELS = [
    "mcsym",
    "vfold",
    "rnacomposer",
    "simrna",
    "3drna",
    "isrna",
    "rhofold",
    "trrosettarna",
    "vfoldpipeline",
    "rnajp",
]
OLD_TO_NEW = {
    "BARNABA-eRMSD": "εRMSD",
    "rnacomposer": "RNAComposer",
    "isrna": "IsRNA1",
    "3drna": "3dRNA",
    "rhofold": "RhoFold",
    "simrna": "SimRNA",
    "vfold": "Vfold3D",
    "rp14_free": "rp14f",
    "rp14_bound": "rp14b",
    "lddt": "lDDT",
    "trrosettarna": "trRosettaRNA",
    "mcsym": "MC-Sym",
    "vfoldpipeline": "Vfold-Pipeline",
    "rnajp": "RNAJP",
}
DESC_METRICS = ["RMSD", "P-VALUE", "εRMSD", "MCQ"]
METRICS = [
    "RMSD",
    "P-VALUE",
    "εRMSD",
    "TM-score",
    "GDT-TS",
    "INF-ALL",
    "CAD",
    "lDDT",
    "MCQ",
]

COLORS = [
    "#e10000",
    "#656567",
    "#ee7f00",
    "#8b1b58",
    "#76885B",
    "#31b2cb",
    "#FFD23F",
    "#B67352",
    "#005793",
]


class PolarChart:
    def __init__(self, in_path: str):
        self.df = self.read_df(in_path, ALL_MODELS, METRICS)

    def _clean_polar_viz(self, fig):
        new_polars = {
            "polar": dict(
                radialaxis=dict(
                    showline=False,
                    showgrid=True,
                    linewidth=0.5,
                    linecolor="black",
                    gridcolor="black",
                    gridwidth=0.5,
                    showticklabels=False,
                ),
                angularaxis=dict(
                    linewidth=0.5,
                    visible=True,
                    linecolor="black",
                    showline=True,
                    gridcolor="black",
                ),
                radialaxis_tickfont_size=20,
                bgcolor="white",
            )
        }
        fig.update_layout(
            legend=dict(
                orientation="v",
                bgcolor="white",
                bordercolor="Black",
                borderwidth=1,
                font=dict(size=20),
                x=1.2,
                y=1.1,
            ),
        )
        fig.update_layout(margin=dict(l=200, r=200, b=50, t=50))
        fig.update_layout(font_size=28)
        fig.update_layout(
            **new_polars,
            showlegend=True,
        )
        return fig

    def viz_dataset(self):
        """Plot the polar distribution for a dataset."""
        df = (
            self.df[["Metric", "Model", "Metric (value)"]]
            .groupby(["Model", "Metric"])
            .mean()
            .reset_index()
        )
        fig = px.bar_polar(
            df,
            r="Metric (value)",
            theta="Model",
            color="Metric",
            template="plotly_white",
            color_discrete_sequence=COLORS,
            range_r=[0, 9],
        )
        fig = self._clean_polar_viz(fig)
        fig.show()
        save_path = "polar_plot.png"
        fig.write_image(save_path, scale=2, width=1000, height=800)

    def normalize_metrics(
        self, df: pd.DataFrame, desc_metrics: List = DESC_METRICS
    ) -> pd.DataFrame:
        """
        Normalize the metrics with the min-max scaler
        :param df: dataframe with the values for each metric
        :param desc_metrics: metrics that are descendant: the lower, the better
        :return: a dataframe with the normalized metrics
        """
        metrics = df["Metric"].unique()
        for metric in metrics:
            mask = df["Metric"] == metric
            metric_values = df.loc[mask, "Metric (value)"].values.reshape(-1, 1)
            if metric_values.shape[0] == 0:
                continue
            non_nan_metrics = metric_values[~np.isnan(metric_values)].reshape(-1, 1)
            scaler = MinMaxScaler().fit(X=non_nan_metrics)
            norm_metric = scaler.transform(X=metric_values).reshape(-1).tolist()
            if metric in desc_metrics:
                norm_metric = [x if np.isnan(x) else 1 - x for x in norm_metric]
            df.loc[mask, "Metric (value)"] = norm_metric
        return df

    def read_df(self, in_path: str, models: List, metrics: List) -> pd.DataFrame:
        """
        Read the dataset results and merge it into one dataframe
        :param in_path: path to a directory with the metric results: one file per RNA
        :param models: list of predictive models to consider
        :param metrics: list of metrics to consider
        :return a dataframe with three column: the metric, the metric name and the model
        """
        out_df = {"Metric": [], "Metric (value)": [], "Model": []}
        # List all the .csv files with the results per RNA
        files = [x for x in os.listdir(in_path) if x.endswith(".csv")]
        for n_file in files:
            in_df = os.path.join(in_path, n_file)
            df = pd.read_csv(in_df, index_col=[0])
            df = df.rename(columns=OLD_TO_NEW)
            # Loop over each predictive model
            for name, row in df.iterrows():
                # Get the name of the given model
                model = name.split("_")[1]
                if model in models:
                    out_df["Metric (value)"].extend(row[metrics])
                    out_df["Metric"].extend(metrics)
                    out_df["Model"].extend([model] * len(metrics))
        out_df = pd.DataFrame(out_df)
        # Function that normalize each metric by the min max scaler
        out_df = self.normalize_metrics(out_df)
        out_df = out_df.replace(OLD_TO_NEW)
        return out_df


if __name__ == "__main__":
    in_path = os.path.join("data", "output")
    polar_chart = PolarChart(in_path)
    polar_chart.viz_dataset()
