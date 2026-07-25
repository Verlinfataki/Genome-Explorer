from genes import rechercher_gene
from collections import Counter

# ----------------------------------------------------------
def traduire_gene(genome_record, nom_gene):
    """
    Traduire un gène en protéine.

    Paramètres
    ----------
    genome_record : SeqRecord
        Génome annoté contenant tous les gènes.

    nom_gene : str
        Nom du gène à traduire.

    Retour
    ------
    bool
        True si le gène existe.
        False sinon.

    Explication biologique
    ----------------------
    Cette fonction extrait la séquence ADN d'un gène
    puis la traduit en une séquence d'acides aminés.

    La traduction suit le code génétique standard.
    """
    
    # ------------------------------------------------------
    # Rechercher le gène dans le génome.
    # Cette fonction parcourt le génome et retourne
    # directement le gène correspondant au nom fourni.
    # ------------------------------------------------------
    gene_feature = rechercher_gene(genome_record, nom_gene)

    # ------------------------------------------------------
    # Vérifier que le gène existe.
    # Si aucun gène n'a été trouvé, la fonction
    # s'arrête immédiatement.
    # ------------------------------------------------------
    if gene_feature is None:

        print(f"\nLe gène '{nom_gene}' est introuvable.")

        return None
        
    # --------------------------------------------------
    # Extraction automatique de la séquence ADN.
    # IMPORTANT :
    # Biopython tient automatiquement compte
    # du brin + ou du brin -.
    # --------------------------------------------------
    sequence_adn = gene_feature.extract(genome_record.seq)

    # --------------------------------------------------
    # Traduction ADN → protéine.
    # to_stop=True signifie :
    # arrêter la traduction dès le premier codon STOP.
    # --------------------------------------------------
    proteine = sequence_adn.translate(to_stop=True)

    # --------------------------------------------------
    # Affichage.
    # --------------------------------------------------


    # Retourner la protéine afin qu'elle puisse être
    # utilisée par d'autres fonctions du projet.
    return proteine

   


# ----------------------------------------------------------
def analyser_composition_proteine(proteine):
    """
    Analyse la composition en acides aminés d'une protéine.

    Paramètre
    ---------
    proteine : Seq ou str
        Séquence protéique à analyser.

    Retour
    ------
    Counter
        Dictionnaire contenant le nombre d'occurrences
        de chaque acide aminé.

    Exemple
    --------
    Entrée :
        MKRISTTITTTITITTGNGAG

    Sortie :
        {'M': 1, 'K': 1, 'R': 1, ...}
    """

    # Transformer la séquence en chaîne de caractères.
    # Cette étape garantit que Counter pourra parcourir
    # correctement chaque acide aminé.
    proteine = str(proteine)

    # Compter automatiquement le nombre d'occurrences
    # de chaque acide aminé.
    composition = Counter(proteine)

    # Retourner le résultat sans l'afficher.
    return composition





