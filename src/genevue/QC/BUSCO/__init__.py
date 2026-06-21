"""
BUSCO - Benchmarking Universal Single-Copy Orthologs.

------------------------------------------------------------------------------------------------------------------------
This is a port of BUSCO v6.0 to GeneVue.
This port was undertaken to facilitate integration with other GeneVue components,
and to avoid over-reliance on the Conda ecosystem.

Major changes:

- All external program calls now go through GeneVue’s custom backend.
- Switched to GeneVue’s log manager.
- Made adjustments to comply with static analysis standards while minimizing changes to program behavior.
- Unlike the original BUSCO, this version does not include any default datasets to reduce the package size. You can
  download the datasets using the command `genevue qc BUSCO --download`.

If you are using only this applet, you need only cite the BUSCO tool and its dataset, not GeneVue.
However, if you are using this applet in conjunction with other components of GeneVue (code outside the `BUSCO` module
but inside GeneVue package), you must cite both the BUSCO tool and its dataset, as well as GeneVue.

Please cite after you using BUSCO and OrthoDB.

For BUSCO v6.0 and OnthoDB v12:

Tegenfeldt F., Kuznetsov D., Manni M., Berkeley M., Zdobnov E.M., Kriventseva E.V.,
OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes.
Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D516–D522, https://doi.org/10.1093/nar/gkae987

For OnthoDB v10:

Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov,
BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic,
Prokaryotic, and Viral Genomes.
Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654, https://doi.org/10.1093/molbev/msab199

BUSCO is under the MIT License. The original license text is located in the `LICENCE_INTEGRATION` file in the root of the
GeneVue repository.

------------------------------------------------------------------------------------------------------------------------
== A Simple User Guide ==

The `BUSCO` applet works exactly like the original BUSCO, but includes some APIs to facilitate interaction with other
GeneVue modules.

Base Command

```
genevue qc BUSCO -i [SEQUENCE_FILE] -m [MODE] [OTHER_OPTIONS]
```

For a description of `-m` (`--mode`), see the "Mode" section; And for `-i` (`--in`), see the “Input or Batch Input” section.

For other parameters:

- `-l` (`--lineage_dataset`), the name of the BUSCO lineage dataset to be used. If you are using a local database, please
  provide its absolute path.
- `-c` (`--cpu`), the number of threads/cores to use, optional. Only using 1 core by default.
- `-o` (`--out`), a recognizable short name for your analyzing run, optional. Using `BUSCO_<input_filename>` by default.

To avoid the hassle of repeatedly entering parameters at the command line, you can provide a general configuration file
to BUSCO. Use the `--config` option to specify the configuration file:

```
genevue qc BUSCO --config /path/to/myconfig.ini
```

Or add an environment variable named `BUSCO_CONFIG_FILE`. For example:

```
export BUSCO_CONFIG_FILE="/path/to/myconfig.ini"
```

Alternatively, enable or adjust the default BUSCO configuration file built into GeneVue (still in .ini format. Recommend):

```
genevue config BUSCO activate    # Make BUSCO use the configuration provided by GeneVue ("configure hijacking")
genevue config BUSCO edit        # Manually edit the configuration file
genevue config BUSCO deactivate  # Prevent GeneVue from "hijacking" BUSCO configurations
```

Modes

BUSCO supports the following these modes (`--mode` or `-m`):
- Genome mode (`--mode genome`), suitable for scenarios where only the genome sequence is available. A gene predictor is
  required (the default is Metaeuk, but Miniport, Augustus, or Prodigal can be used as needed).
- Transcriptome mode (`--mode transcriptome`), which requires the provision of CDS or RNA FASTA file.
- Protein mode (`--mode proteins`), which requires the provision of protein FASTA file.

Structure of the Results Folder



Interpretation of Results

The results will be displayed as a string or, more specifically, in a table:

```
# BUSCO version is: 5.2.2
# The lineage dataset is: saccharomycetes_odb10 (Creation date: 2020-08-05, number of genomes: 76, number of BUSCOs: 2137)
# Summarized benchmarking in BUSCO notation for file /data/manni/busco_protocol/protocol1/Tglobosa_GCF_014133895.1_genome.fna
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk
    Results:
    C:99.6%[S:99.5%,D:0.1%],F:0.1%,M:0.3%,n:2137
    2129    Complete BUSCOs (C)
    2126    Complete and single-copy BUSCOs (S)
    3       Complete and duplicated BUSCOs (D)
    3       Fragmented BUSCOs (F)
    5       Missing BUSCOs (M)
    2137    Total BUSCO groups searched
Dependencies and versions:
    hmmsearch: 3.1
    metaeuk: 4.a0f584d
```

The result consist of two parts:

- The OrthoDB which used (in this case, `saccharomycetes_odb10`) and its release information (`Creation date: 2020-08-05,
  number of genomes: 76, number of BUSCOs: 2137`). The closer the lineage of the database is to the lineage of the genomes,
  the more reliable the results.
- A string containing percentages (`C:99.6%[S:99.5%,D:0.1%],F:0.1%,M:0.3%,n:2137`), the meaning of which is described
  below this string. For a well-assembled genome, `C` should be at least 90% (or even 95%), And `D` should be as low as
  possible.

Plot

Due to compatibility issues, v6 no longer provides the plotting feature (`--plot`); however, as a transitional measure,
GeneVue has written plotting code for that, so the `--plot` feature remains available. When plotting, the program will
still directly read the `short_summary.*.json` file(s) in the specified directory to generate horizontal stacked histograms.

Please note that, due to differences in implementation, the figure formatting may differ from that of v5.
"""

__version__ = "6.0.0"

__all__ = [
    "Actions",
    "Analysis",
    "AutoLineage",
    "augustus",
    "base",
    "bbtools",
    "blast",
    "BuscoAnalysis",
    "BuscoConfig",
    "BuscoDownloadManager",
    "BuscoLogger",
    "BuscoPlacer",
    "BuscoRunner",
    "ConfigManager",
    "Exceptions",
    "GeneSetAnalysis",
    "GenomeAnalysis",
    "hmmer",
    "metaeuk",
    "miniprot",
    "prodigal",
    "run_BUSCO",
    "sepp",
    "Toolset",
    "TranscriptomeAnalysis",
    "__version__",
]
