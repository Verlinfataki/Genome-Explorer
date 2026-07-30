"""
=========================================================
Genome Explorer

Protein analysis
=========================================================
"""

from proteins import (
    calculer_masse,
    calculer_pi,
    calculer_hydrophobicite,
    traduire_gene
)



def analyser_proteines(genome_record):
    """
    Analyse all protein-coding genes.

    Returns
    -------
    list[dict]
    """

    proteins = []

    for feature in genome_record.features:

        if feature.type != "CDS":
            continue



        gene = feature.qualifiers.get("gene", [""])[0]

        if not gene:
            continue

        protein = traduire_gene(genome_record, gene)

        if protein is None:
            continue

        proteins.append(
            {
                "Gene": gene,
                "Protein Length": len(protein),
                "Molecular Weight": calculer_masse(str(protein)),
                "Isoelectric Point": calculer_pi(str(protein)),
                "Hydrophobicity": calculer_hydrophobicite(str(protein)),
            }
        )

    return proteins