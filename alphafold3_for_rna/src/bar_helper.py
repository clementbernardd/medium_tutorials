import os
from typing import Dict, List, Any
import pandas as pd
import numpy as np

import plotly.express as px

from sklearn.preprocessing import MinMaxScaler

CASP_FILES = [
    f"{x}.csv"
    for x in [
        "R1107",
        "R1108",
        "R1116",
        "R1117",
        "R1149",
        "R1156",
        "R1189",
        "R1190",
    ]
]
SUB_METRICS = [
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
    "alphafold3",
    "best",
]

OLD_TO_NEW = {
    "BARNABA-eRMSD": "εRMSD",
    "rnacomposer": "RNAComposer (TP)",
    "isrna": "IsRNA1 (AI)",
    "3drna": "3dRNA (TP)",
    "rhofold": "RhoFold (DL)",
    "simrna": "SimRNA (AI)",
    "vfold": "Vfold3D (TP)",
    "rp14_free": "rp14f",
    "rp14_bound": "rp14b",
    "lddt": "lDDT",
    "trrosettarna": "trRosettaRNA (DL)",
    "mcsym": "MC-Sym (TP)",
    "vfoldpipeline": "Vfold-Pipeline (TP)",
    "rnajp": "RNAJP (AI)",
    "alphafold": "AlphaFold 3",
    "alphafold3": "AlphaFold 3",
    "best": "Challenge-best",
    "RNA_PUZZLES": "RNA-Puzzles",
    "RNASOLO": "RNASolo",
    "RNA3DB_LONG": "RNA3DB_Long",
    "RNA3DB": "RNA3DB_0",
    "CASP_RNA": "CASP-RNA",
    "CASP": "CASP-RNA",
}
COLORS_MAPPING = {
    "RMSD": "#e10000",
    "INF-ALL": "#656567",
    "CAD": "#ee7f00",
    "TM-score": "#8b1b58",
    "GDT-TS": "#76885B",
    "lDDT": "#31b2cb",
    "P-VALUE": "#B67352",
    "εRMSD": "#FFD23F",
    "MCQ": "#005793",
    "INF-STACK": "#ef8927",
    "INF-WC": "#83b8d6",
    "INF-NWC": "#621038",
}
DESC_METRICS = ["RMSD", "P-VALUE", "DI", "εRMSD", "MCQ"]
ORDER_MODELS = ["MC-Sym", "Vfold3D", "RNAComposer", "3dRNA", "SimRNA", "IsRNA1",
                "RhoFold", "trRosettaRNA", "Vfold-Pipeline", "RNAJP", "AlphaFold3", "Challenge-best"]


def update_bar_plot(fig: Any) -> Any:
    """
    Update the bar plot to clean it up
    :param fig: the plotly figure
    """
    params_axes = dict(
        showgrid=True,
        gridcolor="#d6d6d6",
        linecolor="black",
        zeroline=False,
        linewidth=1,
        showline=True,
        mirror=True,
        gridwidth=1,
        griddash="dot",
    )
    fig.update_xaxes(**params_axes)
    fig.update_yaxes(**params_axes)
    fig.update_layout(dict(plot_bgcolor="white"), margin=dict(l=0, r=5, b=0, t=20))
    fig.update_layout(
        font=dict(
            family="Computer Modern",
            size=26,
        )
    )
    fig.update_layout(
        legend=dict(
            orientation="h",
            bgcolor="#f3f3f3",
            bordercolor="black",
            borderwidth=1,
            x=-0.12,
            y=-0.25,
        ),
    )
    return fig


class BarHelper:
    def __init__(self, in_paths: Dict):
        """
        Class that does the bar plot visualisation
        :param in_paths: dictionary with the paths to the .csv files for each dataset
        """
        self.df = self.read_df(in_paths)

    def get_mean_metrics(
        self, in_path: str, metrics: List = SUB_METRICS, models: List = ALL_MODELS
    ) -> Dict:
        """
        Return the mean per metric from a directory with .csv files
        :param in_path: path to a directory with .csv files
        :param metrics: list of metrics to consider
        :param models: list of models to consider
        :return: a dictionary with, as keys, the models and as values the metrics values
        """
        files = [x for x in os.listdir(in_path) if x.endswith(".csv")]
        scores = {model: {metric: [] for metric in metrics} for model in models}
        # Compare CASP-RNA with only structures where almost all the models made a prediction
        if "casp" in in_path:
            files = CASP_FILES
            models = models.copy()
            # Remove MC-Sym from the models
            models.remove("mcsym")
        for n_file in files:
            in_df = os.path.join(in_path, n_file)
            df = pd.read_csv(in_df, index_col=[0])
            for model in models:
                # Get the metrics for the model
                out = self.get_metrics_from_model(df, model, metrics)
                for c_score, metric in zip(out, metrics):
                    scores[model][metric].append(c_score)
        # Get the mean values per benchmark
        scores = self.get_mean_scores(scores)
        return scores

    def get_mean_scores(self, scores: Dict) -> Dict:
        """
        Return the mean scores for each model and metric per benchmark
        :param scores: dictionary with the metrics for each model
        :return: a dictionary with the mean metrics for each model
        """
        for model, values in scores.items():
            for metric_name, metric in values.items():
                scores[model][metric_name] = np.nanmean(metric)
        return scores

    def read_df(self, in_paths: Dict) -> pd.DataFrame:
        """
        Read the dataset results and store it into one common dataframe.
        :param in_paths: dictionary with the paths to the .csv files for each dataset
        :return a dataframe with the metrics for each model and each dataset
        """
        df = {"Metric": [], "Dataset": [], "Metric (value)": [], "Model": []}
        # Loop over each dataset
        for dataset, d_path in in_paths.items():
            # Get the mean values for each model for this given
            c_scores = self.get_mean_metrics(d_path)
            for model, values in c_scores.items():
                metric = list(values.values())
                n = len(metric)
                metric_name = list(values.keys())
                df["Metric"].extend(metric_name)
                df["Dataset"].extend([dataset] * n)
                df["Metric (value)"].extend(metric)
                df["Model"].extend([model] * n)
        df = pd.DataFrame(df)
        # Convert the metric to min-max scale
        df = self.normalize_metrics(df)
        df = df.replace(OLD_TO_NEW)
        return df

    def normalize_metrics(self, df, desc_metrics: List = DESC_METRICS) -> pd.DataFrame:
        """
        Normalize the metrics to a min-max scale
        :param df: dataframe with the metrics for each model and each dataset
        :param desc_metrics: the descending metrics (the lower, the better)
        :return: a dataframe with the normalized metrics
        """
        metrics, datasets = df["Metric"].unique(), df["Dataset"].unique()
        for metric in metrics:
            mask = df["Metric"] == metric
            metric_values = df.loc[mask, "Metric (value)"].values.reshape(-1, 1)
            if metric_values.shape[0] == 0:
                continue
            non_nan_metrics = metric_values[~np.isnan(metric_values)].reshape(-1, 1)
            if len(non_nan_metrics) == 0:
                continue
            scaler = MinMaxScaler().fit(X=non_nan_metrics)
            norm_metric = scaler.transform(X=metric_values).reshape(-1).tolist()
            if metric in desc_metrics:
                norm_metric = [x if np.isnan(x) else 1 - x for x in norm_metric]
            df.loc[mask, "Metric (value)"] = norm_metric
        return df

    def viz(self):
        datasets = self.df["Dataset"].unique()
        for i, dataset in enumerate(datasets):
            self.viz_dataset(dataset)

    def viz_dataset(self, dataset: str):
        """Plot the polar distribution for a dataset."""
        width, height = 1000, 600
        df = (
            self.df[["Metric", "Model", "Metric (value)"]]
            .groupby(["Model", "Metric"])
            .mean()
            .reset_index()
        )
        fig = px.bar(
            df,
            y="Model",
            x="Metric (value)",
            color="Metric",
            color_discrete_map=COLORS_MAPPING,
            orientation="h",
            category_orders={
                "Metric": list(COLORS_MAPPING.keys()),
                "Model": ORDER_MODELS,
            },
            labels={"Metric (value)": "Normalized metrics"},
            range_x=[0, 9],
        )
        fig = update_bar_plot(fig)
        fig.show()
        save_path = os.path.join("data", "plots", dataset + ".png")
        fig.write_image(save_path, scale=2, width=width, height=height)

    def get_metrics_from_model(
        self, df: pd.DataFrame, model: str, metrics: List
    ) -> List:
        """
        Return the metrics for a givel model
        :param df: a dataframe with the metrics for a given RNA
        :param model: the model to return the metric values
        :param metrics: the given metrics to consider
        :return: a list of metrics for the given model
        """
        names = [x for x in df.index if model in x]
        df.rename(columns=OLD_TO_NEW, inplace=True)
        df = df.loc[names].mean(axis=0)
        output = []
        for metric in metrics:
            if metric in df:
                output.append(df[metric])
            else:
                output.append(np.nan)
        return output


if __name__ == "__main__":
    benchmarks = ["CASP_RNA", "RNA_PUZZLES", "RNASOLO", "RNA3DB", "RNA3DB_LONG"]
    in_paths = {name: os.path.join("data", "output", name) for name in benchmarks}
    viz_polar = BarHelper(in_paths)
    viz_polar.viz()
