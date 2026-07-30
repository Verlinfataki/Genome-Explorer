"""
=========================================================
Genome Explorer

Scientific Figure Generator
=========================================================
"""
from fileinput import filename
from pathlib import Path
import matplotlib.pyplot as plt



class PlotGenerator:
    """
    Generate scientific figures.
    """

    def __init__(self, output_dir="outputs/figures"):

        self.output_dir = Path(output_dir)

        self.output_dir.mkdir(
            parents=True,
            exist_ok=True
        )

    def save_figure(self, filename):
        """
        Save current figure.
        """

        output_file = self.output_dir / filename

        plt.tight_layout()

        plt.savefig(
            output_file,
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

        print(f"[OK] {output_file}")


    def histogram(
        self,
        dataframe,
        column,
        title,
        xlabel,
        filename,
        bins=30,
    ):
        """
        Generate a histogram from a dataframe column.
        """

        plt.figure(figsize=(8, 5))

        plt.hist(
            dataframe[column],
            bins=bins,
        )

        plt.title(title)

        plt.xlabel(xlabel)

        plt.ylabel("Frequency")

        self.save_figure(filename)


    def gene_length_distribution(self, genes_dataframe):
            """
            Plot the distribution of gene lengths.
            """
    
            self.histogram(
                genes_dataframe,
                column="Length",
                title="Gene Length Distribution",
                xlabel="Length (bp)",
                filename="gene_length_distribution.png",
            )


    def protein_length_distribution(self, proteins_dataframe):
        """
        Plot protein length distribution.
        """

        self.histogram(
            proteins_dataframe,
            column="Protein Length",
            title="Protein Length Distribution",
            xlabel="Length (amino acids)",
            filename="protein_length_distribution.png",
        )


    def molecular_weight_distribution(self, proteins_dataframe):
        """
        Plot molecular weight distribution.
        """

        self.histogram(
            proteins_dataframe,
            column="Molecular Weight",
            title="Molecular Weight Distribution",
            xlabel="Molecular Weight (Da)",
            filename="molecular_weight_distribution.png",
        )


    def isoelectric_point_distribution(self, proteins_dataframe):
        """
        Plot isoelectric point distribution.
        """

        self.histogram(
            proteins_dataframe,
            column="Isoelectric Point",
            title="Isoelectric Point Distribution",
            xlabel="Isoelectric Point (pI)",
            filename="isoelectric_point_distribution.png",
        )


    def hydrophobicity_distribution(self, proteins_dataframe):
        """
        Plot hydrophobicity distribution.
        """

        self.histogram(
            proteins_dataframe,
            column="Hydrophobicity",
            title="Hydrophobicity Distribution",
            xlabel="Hydrophobicity (%)",
            filename="hydrophobicity_distribution.png",
        )


    def bar_chart(
        self,
        labels,
        values,
        title,
        xlabel,
        ylabel,
        filename,
    ):
        """
        Generate a bar chart.
        """

        plt.figure(figsize=(6, 5))

        plt.bar(labels, values)

        plt.title(title)

        plt.xlabel(xlabel)

        plt.ylabel(ylabel)

        self.save_figure(filename)


    def strand_distribution(self, genes_dataframe):
        """
        Plot strand distribution.
        """

        counts = genes_dataframe["Strand"].value_counts()

        self.bar_chart(
            labels=counts.index.astype(str),
            values=counts.values,
            title="Gene Strand Distribution",
            xlabel="DNA Strand",
            ylabel="Number of Genes",
            filename="strand_distribution.png",
        )