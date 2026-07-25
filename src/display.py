from genes import rechercher_gene

# AFFICHAGE DES RÉSULTATS
def afficher_resume(genome_record, compteur):
    """
    Afficher un résumé des annotations du génome.
    Paramètres
    ----------
    genome_record : SeqRecord
    compteur : Counter
    Retour
    ------
    Aucun.
    """

    print("=" * 60)
    print("EXPLORATION DES GÈNES")
    print("=" * 60)

    print(f"Accession : {genome_record.id}")
    print(f"Description : {genome_record.description}")


    print("-" * 60)

    print(f"Nombre total d'annotations : {len(genome_record.features)}")

    print("-" * 60)

    print("Types d'annotations présents :")

    # Affichage trié par fréquence décroissante.
    for feature_type, nombre in compteur.most_common():

        print(f"{feature_type:<20} : {nombre}")

    print("=" * 60)


# Affichage des 10 premiers gènes
def afficher_premiers_genes(genome_record, limite=10):
    """
    Afficher les premiers gènes du génome.
    Paramètres
    ----------
    genome_record : SeqRecord
        Génome annoté.
    limite : int
        Nombre maximum de gènes à afficher.
    Retour
    ------
    Aucun.

    Explication biologique
    ----------------------
    Chaque annotation de type 'gene' représente un gène
    identifié sur le chromosome. Cette fonction affiche
    les premiers gènes afin de vérifier que les annotations
    ont été correctement lues.
    """

    print("\n" + "=" * 60)
    print(f"PREMIERS {limite} GÈNES")
    print("=" * 60)

    compteur = 0

    for feature in genome_record.features:

        # On ignore tout ce qui n'est pas un gène.
        if feature.type != "gene":
            continue

        compteur += 1

        # Nom du gène.
        gene_name = feature.qualifiers.get("gene", ["Inconnu"])[0]

        # Coordonnées.
        start = int(feature.location.start)
        end = int(feature.location.end)

        # Longueur du gène.
        longueur = end - start

        # Brin d'ADN.
        if feature.location.strand == 1:
            strand = "+"
        elif feature.location.strand == -1:
            strand = "-"
        else:
            strand = "?"

        print(f"Gène {compteur}")
        print(f"Nom      : {gene_name}")
        print(f"Début    : {start}")
        print(f"Fin       : {end}")
        print(f"Longueur : {longueur} pb")
        print(f"Brin      : {strand}")
        print("-" * 60)

        # Arrêt lorsque la limite est atteinte.
        if compteur >= limite:
            break



# ----------------------------------------------------------
def afficher_sequence_gene(genome_record, nom_gene):
    """
    Afficher la séquence ADN d'un gène.
    """

    # ------------------------------------------------------
    # Rechercher le gène dans le génome.
    # ------------------------------------------------------
    feature = rechercher_gene(genome_record, nom_gene)

    # ------------------------------------------------------
    # Vérifier que le gène existe.
    # ------------------------------------------------------
    if feature is None:

        print(f"\nLe gène '{nom_gene}' est introuvable.")

        return False

    # ------------------------------------------------------
    # Extraire automatiquement la séquence ADN.
    #
    # Biopython tient compte du brin (+ ou -)
    # sans intervention de notre part.
    # ------------------------------------------------------
    sequence = feature.extract(genome_record.seq)

    # ------------------------------------------------------
    # Affichage.
    # ------------------------------------------------------
    print("\n" + "=" * 60)
    print("SÉQUENCE DU GÈNE")
    print("=" * 60)

    print(f"Nom : {nom_gene}")

    print(f"Longueur : {len(sequence)} pb")

    print()

    print(sequence)

    print("=" * 60)

    return True


# ----------------------------------------------------------
def afficher_statistiques(statistiques):
    """
    Afficher les statistiques descriptives des gènes.

    Paramètres
    ----------
    statistiques : dict
        Dictionnaire retourné par calculer_statistiques_genes().
    """

    print("\n" + "=" * 60)
    print("STATISTIQUES DES GÈNES")
    print("=" * 60)

    print(f"Nombre de gènes      : {statistiques['nombre_genes']}")
    print(f"Longueur minimale    : {statistiques['longueur_min']} pb")
    print(f"Longueur maximale    : {statistiques['longueur_max']} pb")
    print(f"Longueur moyenne     : {statistiques['longueur_moyenne']:.2f} pb")
    print(f"Longueur médiane     : {statistiques['longueur_mediane']:.2f} pb")
    print(f"Gènes sur le brin +  : {statistiques['genes_plus']}")
    print(f"Gènes sur le brin -  : {statistiques['genes_moins']}")

    print("=" * 60)


def afficher_traduction(proteine, nom_gene, sequence_adn):
    """
    Afficher la séquence protéique traduite à partir d'un gène.

    Paramètres
    ----------
    proteine : Seq ou str
        Séquence protéique à afficher.
    nom_gene : str
        Nom du gène traduit.
    sequence_adn : Seq
        Séquence ADN du gène.
    """

    print("\n" + "=" * 60)
    print("TRADUCTION DU GÈNE")
    print("=" * 60)
    
    print(f"Gène            : {nom_gene}")
    
    print(f"Longueur ADN    : {len(sequence_adn)} nucléotides")
    
    print(f"Longueur protéine : {len(proteine)} acides aminés")
    
    print()
    
    print("Séquence protéique :")
    
    print(proteine)
    
    print("=" * 60)


def afficher_composition_proteine(composition):
    """
    Affiche la composition en acides aminés d'une protéine.

    Paramètre
    ---------
    composition : Counter
        Dictionnaire contenant le nombre d'occurrences
        de chaque acide aminé.

    Explication biologique
    ----------------------
    Chaque lettre représente un acide aminé.

    Cette fonction affiche combien de fois
    chaque acide aminé apparaît dans la protéine.
    """

    # Afficher un titre pour rendre les résultats lisibles.
    print("\n============================================================")
    print("COMPOSITION DE LA PROTÉINE")
    print("============================================================")

    # Calculer la longueur totale de la protéine.
    # Elle correspond à la somme des occurrences
    # de tous les acides aminés.
    longueur = sum(composition.values())

    print(f"Longueur : {longueur} acides aminés")
    print()

    # Parcourir les acides aminés par ordre alphabétique
    # afin d'obtenir un affichage toujours identique.
    for acide_amine in sorted(composition):

        print(f"{acide_amine} : {composition[acide_amine]}")