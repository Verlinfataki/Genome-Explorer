from Bio.SeqUtils.ProtParam import ProteinAnalysis
from genes import rechercher_gene
from collections import Counter
from Bio import Align
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

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


    reste = len(sequence_adn) % 3

    if reste != 0:
        sequence_adn = sequence_adn[:-reste]

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



def comparer_traductions(traduction_calculee, traduction_genbank):
    """
    Compare la traduction calculée avec celle fournie par GenBank.

    Paramètres
    ----------
    traduction_calculee : str
        Protéine obtenue avec notre algorithme.
    traduction_genbank : str
        Protéine officielle provenant de GenBank.

    Retour
    ------
    bool
        True si les deux séquences sont identiques,
        False sinon.
    """

    return traduction_calculee == traduction_genbank


def calculer_masse(proteine):
    """
    Calcule la masse moléculaire d'une protéine.

    Paramètre
    ---------
    proteine : str
        Séquence protéique.

    Retour
    ------
    float
        Masse moléculaire en Daltons (Da).
    """

    # Analyse de la protéine
    analyse = ProteinAnalysis(proteine)

    # Retourne la masse moléculaire
    return analyse.molecular_weight()


def calculer_pi(proteine):
    """
    Calcule le point isoélectrique d'une protéine.

    Paramètre
    ---------
    proteine : str
        Séquence protéique.

    Retour
    ------
    float
        Point isoélectrique.
    """

    # Analyse de la protéine
    analyse = ProteinAnalysis(proteine)

    # Retourne le point isoélectrique
    return analyse.isoelectric_point()


def calculer_hydrophobicite(proteine):
    """
    Calcule le pourcentage d'acides aminés hydrophobes.

    Paramètre
    ---------
    proteine : str
        Séquence protéique.

    Retour
    ------
    float
        Pourcentage d'acides aminés hydrophobes.
    """

    # Acides aminés hydrophobes
    hydrophobes = {"A", "V", "I", "L", "M", "F", "W", "Y", "P"}

    # Nombre d'acides aminés hydrophobes
    nb = sum(1 for aa in proteine if aa in hydrophobes)

    # Calcul du pourcentage
    return (nb / len(proteine)) * 100



def comparer_proteines(proteine1, proteine2):
    """
    Compare deux protéines position par position.

    Paramètres
    ----------
    proteine1 : str
        Première séquence protéique.
    proteine2 : str
        Deuxième séquence protéique.

    Retour
    ------
    dict
        Résultats de la comparaison.
    """

    # Longueur commune
    longueur = min(len(proteine1), len(proteine2))

    # Nombre d'acides aminés identiques
    identiques = sum(
        1
        for aa1, aa2 in zip(proteine1[:longueur], proteine2[:longueur])
        if aa1 == aa2
    )

    # Pourcentage d'identité
    similarite = (identiques / longueur) * 100 if longueur > 0 else 0

    return {
        "longueur1": len(proteine1),
        "longueur2": len(proteine2),
        "identiques": identiques,
        "longueur_comparee": longueur,
        "similarite": similarite
    }


def aligner_proteines(proteine1, proteine2):
    """
    Réalise un alignement global de deux protéines.

    Paramètres
    ----------
    proteine1 : str
    proteine2 : str

    Retour
    ------
    Alignment
        Meilleur alignement trouvé.
    """

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"

    alignements = aligner.align(proteine1, proteine2)

    return alignements[0]


def lancer_blast(proteine):
    """
    Lance une recherche BLASTP sur le serveur du NCBI.

    Paramètre
    ---------
    proteine : str
        Séquence protéique.

    Retour
    ------
    Blast
        Résultat BLAST.
    """

    # Envoi de la séquence au NCBI
    resultat = NCBIWWW.qblast(
        "blastp",
        "nr",
        proteine
    )

    # Lecture du résultat XML
    blast = NCBIXML.read(resultat)

    return blast


def meilleurs_hits(blast, limite=5):
    """
    Extrait les meilleurs résultats BLAST.
    """

    hits = []

    for alignement in blast.alignments[:limite]:

        hsp = alignement.hsps[0]

        hits.append({
            "accession": alignement.accession,
            "description": alignement.hit_def.split(">")[0][:120],
            "score": hsp.score,
            "identite": (
                hsp.identities / hsp.align_length
            ) * 100,
            "longueur": hsp.align_length
        })

    return hits