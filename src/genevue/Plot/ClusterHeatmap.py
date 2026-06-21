""" """

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


class CHeatmapDrawer:
    """
    A simple class to draw a customizable and elegant cluster heatmap.
    ```
    """

    def __init__(self):
        self.data = sns.load_dataset("iris")
        self.fig, self.ax = None, None

    def draw(self):
        row_colors = self.data["species"].map(
            {"setosa": "red", "versicolor": "green", "virginica": "blue"}
        )
        self.data.pop("species")
        sns.clustermap(
            data=self.data,
            cmap="vlag",
            z_score=0,
            linewidths=0,
            xticklabels=True,
            yticklabels=False,
            row_colors=row_colors,
            method="average",
            metric="euclidean",
            figsize=(10, 6),
            cbar_kws={"label": "Value"},
        )
        plt.show()
