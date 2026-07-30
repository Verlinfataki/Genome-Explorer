from statistics import mean, median

# ----------------------------------------------------------
def calculer_statistiques_genes(genome_record):
    """
    Calculer les statistiques descriptives des gènes.

    Paramètres
    ----------
    genome_record : SeqRecord
        Génome annoté.

    Retour
    ------
    dict
        Dictionnaire contenant les principales statistiques
        sur les gènes du génome.

    Explication biologique
    ----------------------
    Les statistiques descriptives permettent d'obtenir une
    vue d'ensemble de l'organisation du génome avant toute
    analyse plus approfondie.
    """

    longueurs = []

    genes_plus = 0
    genes_moins = 0

    for feature in genome_record.features:

        if feature.type != "gene":
            continue

        longueur = int(feature.location.end) - int(feature.location.start)
        longueurs.append(longueur)

        if feature.location.strand == 1:
            genes_plus += 1

        elif feature.location.strand == -1:
            genes_moins += 1

    statistiques = {
        "gene_count": len(longueurs),
        "min_length": min(longueurs),
        "max_length": max(longueurs),
        "mean_length": mean(longueurs),
        "median_length": median(longueurs),
        "plus_strand": genes_plus,
        "minus_strand": genes_moins,
    }

    return statistiques


