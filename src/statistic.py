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
        "nombre_genes": len(longueurs),
        "longueur_min": min(longueurs),
        "longueur_max": max(longueurs),
        "longueur_moyenne": mean(longueurs),
        "longueur_mediane": median(longueurs),
        "genes_plus": genes_plus,
        "genes_moins": genes_moins,
    }

    return statistiques


