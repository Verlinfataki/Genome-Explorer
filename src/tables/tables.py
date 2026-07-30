"""
=========================================================
Genome Explorer

Scientific Table Generator
=========================================================
"""

from pathlib import Path

import pandas as pd


class TableGenerator:
    """
    Generate scientific tables (CSV files).
    """

    def __init__(self, output_dir="outputs/tables"):

        self.output_dir = Path(output_dir)

        self.output_dir.mkdir(
            parents=True,
            exist_ok=True
        )

    def save_dataframe(self, dataframe, filename):
        """
        Save a DataFrame as CSV.
        """

        output_file = self.output_dir / filename

        dataframe.to_csv(
            output_file,
            index=False,
            encoding="utf-8"
        )

        print(f"[OK] {output_file}")


    def create_table(self, data, filename):
        """
        Create and save a CSV table from a list of dictionaries.

        Parameters
        ----------
        data : list[dict]
            Data to convert into a DataFrame.

        filename : str
            Output CSV filename.

        Returns
        -------
        pandas.DataFrame
        """

        dataframe = pd.DataFrame(data)

        self.save_dataframe(dataframe, filename)

        return dataframe


    def gene_statistics(self, statistics):
            """
            Generate the gene statistics table.
            """
    
            data = pd.DataFrame([
                {
                    "Metric": "Gene count",
                    "Value": statistics["gene_count"],
                    "Unit": "genes"
                },
                {
                    "Metric": "Minimum length",
                    "Value": statistics["min_length"],
                    "Unit": "bp"
                },
                {
                    "Metric": "Maximum length",
                    "Value": statistics["max_length"],
                    "Unit": "bp"
                },
                {
                    "Metric": "Mean length",
                    "Value": statistics["mean_length"],
                    "Unit": "bp"
                },
                {
                    "Metric": "Median length",
                    "Value": statistics["median_length"],
                    "Unit": "bp"
                },
                {
                    "Metric": "Genes on + strand",
                    "Value": statistics["plus_strand"],
                    "Unit": "genes"
                },
                {
                    "Metric": "Genes on - strand",
                    "Value": statistics["minus_strand"],
                    "Unit": "genes"
                }
            ] )
    
            return self.create_table(
                data,
                "gene_statistics.csv"
            )


    def genes(self, genome_record):
        """
            Generate a table containing all annotated genes.
        """

        data = []

        for feature in genome_record.features:

            if feature.type != "gene":
             continue

        qualifiers = feature.qualifiers

        gene = qualifiers.get("gene", [""])[0]
        locus_tag = qualifiers.get("locus_tag", [""])[0]

        start = int(feature.location.start) + 1
        end = int(feature.location.end)

        length = end - start + 1

        strand = "+" if feature.location.strand == 1 else "-"

        product = qualifiers.get("product", [""])[0]

        data.append(
            {
                "Gene": gene,
                "Locus Tag": locus_tag,
                "Start": start,
                "End": end,
                "Length": length,
                "Strand": strand,
                "Product": product
            }
        )

        return self.create_table(
            data,
            "genes.csv"
        )


    def proteins(self, proteins):
        """
            Generate a table containing all analyzed proteins.
        """

        return self.create_table(
            proteins,
            "proteins.csv"
        )
    
    
    