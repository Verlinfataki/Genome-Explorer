# SeqIO permet de lire les fichiers biologiques.
from Bio import SeqIO

# Path facilite la manipulation des chemins de fichiers.
from pathlib import Path
from collections import Counter


# ------------------------------------------------------
# LOCALISATION DU FICHIER GENBANK
# Chemin vers le fichier GenBank téléchargé précédemment.
GENBANK_FILE = Path("data") / "raw" / "ecoli_genome.gb"

# RÉCUPÉRATION DES ANNOTATIONS BIOLOGIQUES
# Lecture du fichier GenBank, la fonction retourne un objet SeqRecord.
# Explication biologique:
"""
    Contrairement au format FASTA, un fichier GenBank contient
    la séquence ADN mais aussi toutes les annotations
    (gènes, CDS, ARN, protéines, etc.).
"""
def lire_genome():
    return SeqIO.read(GENBANK_FILE, "genbank")


# COMPTAGE DES TYPES D'ANNOTATIONS
# Le fichier GenBank contient une liste appelée "features".
# Chaque élément de cette liste est un objet SeqFeature
# représentant une annotation biologique.
# Explication biologique
"""
    Une annotation représente un élément identifié sur le
    chromosome : gène, CDS, tRNA, rRNA, origine de réplication,
    élément mobile, etc.
"""
def compter_annotations(genome_record):
    compteur = Counter()
    for feature in genome_record.features:
        compteur[feature.type] += 1

    return compteur