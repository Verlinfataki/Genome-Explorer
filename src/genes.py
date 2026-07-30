# 
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



def rechercher_cds(genome, nom_gene):
    """
    Recherche la CDS associée à un gène.

    Paramètres
    ----------
    genome : SeqRecord
        Génome chargé depuis le fichier GenBank.
    nom_gene : str
        Nom du gène recherché (ex : "thrL").

    Retour
    ------
    SeqFeature
        La CDS correspondante si elle existe.
        None sinon.
    """

    # Parcourt toutes les annotations du génome
    for feature in genome.features:

        # On ne s'intéresse qu'aux CDS
        if feature.type != "CDS":
            continue

        # Vérifie que l'annotation possède un nom de gène
        if "gene" not in feature.qualifiers:
            continue

        # Compare le nom recherché
        if feature.qualifiers["gene"][0] == nom_gene:
            return feature

    # Aucune CDS trouvée
    return None


def obtenir_infos_cds(cds):
    """
    Extrait les principales informations d'une CDS.

    Paramètre
    ---------
    cds : SeqFeature
        Annotation CDS.

    Retour
    ------
    dict
        Dictionnaire contenant les informations de la CDS.
    """

    # Si aucune CDS n'est fournie
    if cds is None:
        return None

    # Construction du dictionnaire de résultats
    infos = {
        "gene": cds.qualifiers.get("gene", ["Inconnu"])[0],
        "produit": cds.qualifiers.get("product", ["Inconnu"])[0],
        "protein_id": cds.qualifiers.get("protein_id", ["Inconnu"])[0],
        "debut": int(cds.location.start),
        "fin": int(cds.location.end),
        "longueur": len(cds.location),
        "brin": "+" if cds.location.strand == 1 else "-",
        # Séquence protéique officielle de GenBank
        "translation": cds.qualifiers.get("translation", [""])[0]
    }

    return infos



def genes_sans_cds(genome):
    """
    Retourne la liste des gènes qui ne possèdent pas de CDS.

    Paramètre
    ---------
    genome : SeqRecord
        Génome chargé depuis un fichier GenBank.

    Retour
    ------
    list
        Liste des noms des gènes sans CDS.
    """

    # Ensemble contenant tous les noms de gènes
    genes = set()

    # Ensemble contenant tous les noms de gènes ayant une CDS
    genes_cds = set()

    # Parcours de toutes les annotations
    for feature in genome.features:

        # Ignore les annotations sans nom de gène
        if "gene" not in feature.qualifiers:
            continue

        # Nom du gène
        nom = feature.qualifiers["gene"][0]

        # Si c'est une annotation "gene"
        if feature.type == "gene":
            genes.add(nom)

        # Si c'est une annotation "CDS"
        elif feature.type == "CDS":
            genes_cds.add(nom)

    # Les gènes sans CDS sont ceux présents dans genes
    # mais absents de genes_cds
    return sorted(genes - genes_cds)



def obtenir_type_gene(genome, nom_gene):
    """
    Retourne tous les types d'annotations associés à un gène.

    Paramètres
    ----------
    genome : SeqRecord
        Génome chargé.
    nom_gene : str
        Nom du gène.

    Retour
    ------
    list
        Liste des types d'annotations.
    """

    # Ensemble pour éviter les doublons
    types = set()

    # Parcours des annotations
    for feature in genome.features:

        # Vérifie que le nom du gène existe
        if "gene" not in feature.qualifiers:
            continue

        # Vérifie le nom
        if feature.qualifiers["gene"][0] != nom_gene:
            continue

        # Ignore l'annotation "gene"
        if feature.type == "gene":
            continue

        # Ajoute le type rencontré
        types.add(feature.type)

    # Retourne les types triés
    return sorted(types)