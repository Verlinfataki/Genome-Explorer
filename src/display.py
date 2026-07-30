from genes import obtenir_type_gene, rechercher_gene

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

    print(f"Nombre de gènes      : {statistiques['gene_count']}")
    print(f"Longueur minimale    : {statistiques['min_length']} pb")
    print(f"Longueur maximale    : {statistiques['max_length']} pb")
    print(f"Longueur moyenne     : {statistiques['mean_length']:.2f} pb")
    print(f"Longueur médiane     : {statistiques['median_length']:.2f} pb")
    print(f"Gènes sur le brin +  : {statistiques['plus_strand']}")
    print(f"Gènes sur le brin -  : {statistiques['minus_strand']}")

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



def afficher_cds(infos):
    """
    Affiche les informations d'une CDS.
    """

    if infos is None:
        print("\nAucune CDS trouvée.")
        return

    print("\n" + "=" * 60)
    print("INFORMATIONS DE LA CDS")
    print("=" * 60)

    print(f"Gène       : {infos['gene']}")
    print(f"Produit    : {infos['produit']}")
    print(f"Protein ID : {infos['protein_id']}")
    print(f"Début      : {infos['debut']}")
    print(f"Fin        : {infos['fin']}")
    print(f"Longueur   : {infos['longueur']} pb")
    print(f"Brin       : {infos['brin']}")
    print(f"\nTraduction officielle :")
    print(infos["translation"])

    print("=" * 60)

    

def afficher_validation_traduction(est_valide):
    """
    Affiche le résultat de la comparaison des traductions.
    """

    print("\n" + "=" * 60)
    print("VALIDATION DE LA TRADUCTION")
    print("=" * 60)

    if est_valide:
        print("✓ Notre traduction est identique à celle de GenBank.")
    else:
        print("✗ Les deux traductions sont différentes.")

    print("=" * 60)


def afficher_genes_sans_cds(genome,liste):
    """
    Affiche les gènes ne possédant pas de CDS.
    """

    print("\n" + "=" * 60)
    print("GÈNES SANS CDS")
    print("=" * 60)

    print(f"Nombre : {len(liste)}\n")

    # Affiche les 20 premiers gènes
    for nom in liste[:20]:
        type_gene = obtenir_type_gene(genome, nom)

        # Si aucun type n'a été trouvé
        if not type_gene:
            types_str = "Inconnu"
        else:
            types_str = ", ".join(type_gene)

        # Affichage
        print(f"{nom:<20} ({types_str})")

    # Indique s'il reste d'autres gènes
    if len(liste) > 20:
        print(f"\n... {len(liste)-20} autres gènes")

    print("=" * 60)



def afficher_rapport_proteine(infos_cds, proteine, composition, masse, pi, hydrophobicite):
    """
    Affiche un rapport complet sur une protéine.
    """

    print("\n" + "=" * 60)
    print("RAPPORT D'ANALYSE DE LA PROTÉINE")
    print("=" * 60)

    print(f"Gène              : {infos_cds['gene']}")
    print(f"Produit           : {infos_cds['produit']}")
    print(f"Protein ID        : {infos_cds['protein_id']}")
    print(f"Longueur          : {len(proteine)} aa")
    print(f"Masse moléculaire : {masse:.2f} Da")
    print(f"Point isoélectrique : {pi:.2f}")
    print(f"Hydrophobicité    : {hydrophobicite:.2f} %")

    print("\nComposition :")

    for aa in sorted(composition):
        print(f"  {aa} : {composition[aa]}")

    print("=" * 60)


def afficher_comparaison(resultat):
    """
    Affiche le résultat de la comparaison entre deux protéines.
    """

    print("\n" + "=" * 60)
    print("COMPARAISON DE DEUX PROTÉINES")
    print("=" * 60)

    print(f"Longueur protéine 1 : {resultat['longueur1']} aa")
    print(f"Longueur protéine 2 : {resultat['longueur2']} aa")
    print(f"Positions comparées : {resultat['longueur_comparee']}")
    print(f"Acides aminés identiques : {resultat['identiques']}")
    print(f"Similarité : {resultat['similarite']:.2f} %")

    print("=" * 60)


def afficher_alignement(alignement):
    """
    Affiche le meilleur alignement.
    """

    print("\n" + "=" * 60)
    print("ALIGNEMENT DES PROTÉINES")
    print("=" * 60)

    print(f"Score : {alignement.score:.2f}")
    print(f"Longueur séquence 1 : {alignement.shape[0]}")
    print(f"Longueur séquence 2 : {alignement.shape[1]}")

    print("\nAperçu :\n")
    print(str(alignement)[:600])

    print("\n...")

    print("=" * 60)


def afficher_hits_blast(hits):
    """
    Affiche les meilleurs résultats BLAST.
    """

    print("\n" + "=" * 60)
    print("MEILLEURS RÉSULTATS BLAST")
    print("=" * 60)

    for i, hit in enumerate(hits, 1):

        print(f"\nHit {i}")

        print(f"Accession  : {hit['accession']}")
        print(f"Protéine   : {hit['description']}")
        print(f"Score       : {hit['score']}")
        print(f"Identité    : {hit['identite']:.2f} %")
        print(f"Longueur    : {hit['longueur']} aa")

    print("=" * 60)