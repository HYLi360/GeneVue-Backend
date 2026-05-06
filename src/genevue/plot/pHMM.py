import re
from pathlib import Path
from typing import Literal

import logomaker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from genevue import console
from genevue.constants.AminoAcids import AACmap, BLOSUM62


class HMMPreview:
    def __init__(self, phmm_path: str | Path):
        self.phmm_path: Path = (
            Path(phmm_path) if isinstance(phmm_path, str) else phmm_path
        )
        self.color_scheme: Literal["chemistry"] = "chemistry"

        self._parse_phmm()
        console.log(
            f"pHMM file has been loaded. {self.emission.shape[0]} nodes imported."
        )

        self._calc_probs_entropy()
        console.log(
            "Probabities and Infomation Content have been calculated successfully."
        )

    def _parse_phmm(self):
        """
        A simple parser to read Match Emission (-ln P) from phmm file.
        """
        amino_acids = []
        data = []

        with open(self.phmm_path, "r") as f:
            lines = f.readlines()

        start_reading = False

        for i, line in enumerate(lines):
            line = line.strip()

            # 1. fetch amino acids list
            if line.startswith("HMM") and "A" in line and "Y" in line:
                parts = line.split()
                # HMM A C D ... Y (skip first 'HMM')
                amino_acids = parts[1:]
                continue

            # 2. find the start point of data (after COMPO line)
            if line.startswith("COMPO"):
                start_reading = True
                continue

            if start_reading:
                if re.match(r"^\s*\d+\s+\d+\.", line):
                    parts = line.split()
                    # parts[0] = Node ID
                    # parts[1:21] = Score of every a.a.
                    try:
                        scores = [float(x) for x in parts[1:21]]
                        if len(scores) == 20:
                            data.append(scores)
                    except ValueError:
                        continue

                # break if it reads "//" (The model's over)
                if line.startswith("//"):
                    break

        self.amino_acids = amino_acids
        self.emission = np.array(data)

    def _calc_probs_entropy(self):
        # reduction the emission matrix to original
        emission_original = np.exp(-self.emission)

        # normalization to mitigate accuracy errors
        emission_original = emission_original / emission_original.sum(
            axis=1, keepdims=True
        )

        # turn to DataFrame
        self.emission_original = pd.DataFrame(
            data=emission_original,
            columns=self.amino_acids,
        )

        # calculate the Shannon Entropy
        self.entropy = -np.nansum(
            self.emission_original * np.log2(self.emission_original + 1e-10), axis=1
        )

        # calculate the Information Content
        # 1. SeqLogo: IC = Max entropy - entropy
        #    Max entropy = log2(20) or 4.32
        #    REF: Schneider TD, Stephens RM.
        #    Sequence logos: a new way to display consensus sequences.
        #    Nucleic Acids Res. 1990 Oct 25;18(20):6097-100.
        #    doi: 10.1093/nar/18.20.6097. PMID: 2172928; PMCID: PMC332411.
        self.info_content_seq = np.log2(20) - self.entropy

        # 2. HMMLogo: KL_i = sum_j p_ij * log2(p_ij / bg_j)
        #    REF: Eddy, S.R. et,al.
        #    HMMER: biosequence analysis using profile hidden Markov models.
        #    http://hmmer.org
        bgfreq = np.array(
            [BLOSUM62.get(aa, 0) for aa in self.emission_original.columns]
        ).T
        self.info_content_hmm = np.sum(
            self.emission_original.values
            * np.log2(self.emission_original.values / bgfreq),
            axis=1,
        )

    def draw_emission_heatmap(
        self,
        figsize: tuple = (15, 6),
        title: str = "",
        cmap: Literal["Clustal2", "RasMol"] = "Clustal2",
    ):
        aacmap = AACmap(cmap)

        fig, ax = plt.subplots(figsize=(15, 6))

        threshold = 0.2

        for i, aa in enumerate(aacmap.aalabels):
            positions = np.where(self.emission_original.T.values[i] > threshold)[0]
            values = self.emission_original.T.values[i, positions]

            ax.scatter(
                positions,
                [i] * len(positions),
                s=values * 150,
                color=aacmap.color(aa),
                ec="#5F5F5F",
                alpha=0.6,
            )

        ax.set_yticks(range(len(aacmap.aalabels)))
        ax.set_yticklabels(aacmap.aalabels)
        ax.set_xlabel("Model Node Position")
        ax.set_ylabel("Amino Acid")
        ax.set_xlim(0, len(self.info_content_seq))
        ax.grid(True, alpha=0.3)
        ax.set_title(label=title)
        fig.tight_layout()
        plt.show()

        return fig, ax

    def draw_conservation_curse(
        self,
        figsize: tuple = (15, 6),
        title: str = "",
        color: str = "#2c3e50",
        linewidth: float = 1.5,
    ):
        fig, ax = plt.subplots(figsize=figsize)

        # plot the curse
        ax.plot(
            range(1, len(self.info_content_seq) + 1),
            self.info_content_seq,
            color=color,
            linewidth=linewidth,
        )

        # set title and labels
        ax.set_title(label=title, fontsize=14)
        ax.set_xlabel("Model Node Position", fontsize=12)
        ax.set_ylabel("Infomation Content (bits)", fontsize=12)
        ax.set_xlim(0, len(self.info_content_seq))
        ax.grid(True, alpha=0.3)

        # 3 bits reference line
        ax.add_line(
            Line2D(
                xdata=[0, len(self.info_content_seq)],
                ydata=[3, 3],
                linewidth=0.5,
                color="#2c3e50",
            )
        )

        # labeling high-conservation points
        high_cons = np.where(self.info_content_seq > 3.0)[0]
        ax.scatter(
            high_cons + 1,
            self.info_content_seq[high_cons],
            color="red",
            s=20,
            label="High Conservation (> 3 bits)",
        )
        ax.legend()

        fig.tight_layout()

        return fig, ax

    def draw_seqlogo(
        self,
        method: Literal["SeqLogo", "HMMLogo"] = "HMMLogo",
        warpping=50,
        title: str = "",
        highlight=None,
    ):
        console.log("Start drow Sequence Logo.")
        console.log("This process may take some time. Please be patient.")

        # Calculate Infomation Content per amino acid resides
        df = None
        match method:
            case "SeqLogo":
                df = self.emission_original.multiply(self.info_content_seq, axis=0)
            case "HMMLogo":
                df = self.emission_original.multiply(self.info_content_hmm, axis=0)

        # Calculate the seqlogo minimal lines and DF rows
        n_rows = int(np.ceil(len(df) / warpping))
        df_rows = n_rows * warpping

        # Fill 0 to adjust, preventing uneven line lengths
        df = df.reset_index(drop=True)
        df = df.reindex(range(df_rows), fill_value=0.0)
        df.index = np.arange(1, len(df) + 1)

        fig, axes = plt.subplots(n_rows, 1, figsize=(12, 2.5 * n_rows), sharey=True)

        if n_rows == 1:
            axes = [axes]

        max_y = df.sum(axis=1).max() * 1.1

        # 3. 循环绘制每一行
        for i in range(n_rows):
            start = i * warpping
            end = min((i + 1) * warpping, len(df))

            # 切片 DataFrame
            df_slice = df.iloc[start:end]

            # 在当前子图 (axes[i]) 上绘制 Logo
            logo = logomaker.Logo(
                df_slice,
                ax=axes[i],
                color_scheme="chemistry",
                vpad=0.05,
                width=0.8,
            )

            # 样式美化
            logo.style_spines(visible=False)
            logo.style_spines(spines=("left", "bottom"), visible=True)

            # Y
            axes[i].set_ylim(0, max_y)
            axes[i].set_ylabel("Info Content (bits)")

            # 设置 X 轴刻度 (关键：要显示真实的序列位置)
            # logomaker 默认从索引开始，如果你的索引已经是真实位置则无需调整
            # 这里假设 df.index 是 1, 2, 3...
            axes[i].set_xlim(df_slice.index[0] - 0.5, df_slice.index[-1] + 0.5)

            # high-light motif
            if highlight:
                for pos in highlight:
                    if start < pos <= end:
                        # 在背景画一个浅红色竖条
                        axes[i].axvspan(pos - 0.5, pos + 0.5, color="red", alpha=0.2)
                        # 可以在上方加个标记
                        axes[i].text(
                            pos, max_y * 0.9, "▼", ha="center", color="red", fontsize=10
                        )

        # Final adjusting
        fig.suptitle(title, fontsize=16, y=0.99)
        fig.tight_layout()

        return fig, axes
