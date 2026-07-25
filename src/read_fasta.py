"""
Mission 3 : Lire un fichier FASTA

Objectif :
Découvrir le contenu d'un fichier FASTA
grâce à Biopython.

Auteur : Verlin Fataki
Projet : Genome Explorer
"""

# ------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------

# Permet de lire les fichiers biologiques (FASTA, GenBank, etc.).
from Bio import SeqIO

# Permet de manipuler les chemins de fichiers proprement.
from pathlib import Path


# ------------------------------------------------------
# LOCALISATION DU FICHIER FASTA
# ------------------------------------------------------

# Chemin vers le génome téléchargé lors de la mission 2.
FASTA_FILE = Path("data") / "raw" / "ecoli_genome.fasta"

# ------------------------------------------------------
# LECTURE DU FICHIER FASTA
# ------------------------------------------------------

# Lire le fichier FASTA.
# La fonction SeqIO.read() retourne un objet SeqRecord.

genome_record = SeqIO.read(FASTA_FILE, "fasta")


# ------------------------------------------------------
# IDENTIFIANT DU GÉNOME
# ------------------------------------------------------

# L'attribut 'id' contient l'identifiant de la séquence.
print(f"Identifiant : {genome_record.id}")

# ------------------------------------------------------
# DESCRIPTION DU GÉNOME
# ------------------------------------------------------

# La description correspond à toute la première ligne du FASTA.
print(f"Description : {genome_record.description}")

# ------------------------------------------------------
# SÉQUENCE D'ADN
# ------------------------------------------------------

# L'attribut 'seq' contient la séquence biologique.
#print(genome_record.seq)

# ------------------------------------------------------
# LONGUEUR DU GÉNOME
# ------------------------------------------------------

# La fonction len() retourne le nombre de nucléotides.
genome_length = len(genome_record.seq)

print(f"Longueur du génome : {genome_length} nucléotides")