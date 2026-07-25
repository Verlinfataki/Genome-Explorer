from Bio import SeqIO

# ----------------------------------------------------------
def rechercher_gene(genome_record, nom_gene):
    """
    Rechercher un gène dans le génome.

    Paramètres
    ----------
    genome_record : SeqRecord
        Génome annoté.

    nom_gene : str
        Nom du gène recherché.

    Retour
    ------
    SeqFeature | None
        Retourne l'objet représentant le gène.
        Retourne None si le gène n'existe pas.

    Explication biologique
    ----------------------
    Un génome peut contenir plusieurs milliers de gènes.

    Cette fonction joue le rôle d'un moteur de recherche :
    elle parcourt toutes les annotations jusqu'à trouver
    le gène demandé.
    """

    # ------------------------------------------------------
    # Parcourir toutes les annotations du génome.
    # ------------------------------------------------------
    for feature in genome_record.features:

        # Nous nous intéressons uniquement aux annotations
        # de type "gene".
        
        if feature.type != "gene":
            continue

        # --------------------------------------------------
        # Récupérer le nom du gène.
        #
        # Certains gènes peuvent ne pas posséder
        # de champ "gene".
        # --------------------------------------------------
        gene_name = feature.qualifiers.get("gene", [""])[0]

        # --------------------------------------------------
        # Si le nom correspond,
        # on retourne immédiatement ce gène.
        #
        # Pourquoi ?
        #
        # Parce qu'un gène trouvé n'a plus besoin
        # que la boucle continue.
        # --------------------------------------------------
        if gene_name == nom_gene:
            return feature

    # ------------------------------------------------------
    # Aucun gène trouvé.
    # ------------------------------------------------------
    return None



