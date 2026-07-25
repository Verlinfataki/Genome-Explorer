"""
Mission 4 : Exploration du génome

Objectif :
Analyser la composition nucléotidique d'un génome
et calculer les principales statistiques biologiques.

Auteur : Verlin Fataki
Projet : Genome Explorer
"""

# ------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------

from Bio import SeqIO
from pathlib import Path


# ------------------------------------------------------
# LOCALISATION DU FICHIER FASTA

FASTA_FILE = Path("data") / "raw" / "ecoli_genome.fasta"

# ------------------------------------------------------
# LECTURE DU GÉNOME

genome_record = SeqIO.read(FASTA_FILE, "fasta")

# La séquence est convertie en chaîne de caractères
# afin de faciliter les opérations de comptage.
genome_sequence = str(genome_record.seq)

# ------------------------------------------------------
# INFORMATIONS GÉNÉRALES

genome_accession = genome_record.id
genome_description = genome_record.description
genome_length = len(genome_sequence)

# ------------------------------------------------------
# COMPOSITION NUCLÉOTIDIQUE
# ------------------------------------------------------

adenine_count = genome_sequence.count("A")
thymine_count = genome_sequence.count("T")
guanine_count = genome_sequence.count("G")
cytosine_count = genome_sequence.count("C")

# ------------------------------------------------------
# POURCENTAGES
# ------------------------------------------------------

gc_count = guanine_count + cytosine_count
at_count = adenine_count + thymine_count

gc_content = (gc_count / genome_length) * 100
at_content = (at_count / genome_length) * 100

# ------------------------------------------------------
# NUCLÉOTIDE LE PLUS FRÉQUENT
# ------------------------------------------------------

nucleotide_counts = {
    "A": adenine_count,
    "T": thymine_count,
    "G": guanine_count,
    "C": cytosine_count,
}

most_frequent_nucleotide = max(
    nucleotide_counts,
    key=nucleotide_counts.get
)

# ------------------------------------------------------
# AFFICHAGE
# ------------------------------------------------------

print("=" * 60)
print("EXPLORATION DU GÉNOME")
print("=" * 60)

print(f"Accession : {genome_accession}")
print(f"Description : {genome_description}")
print(f"Longueur : {genome_length:,} nucléotides")

print("-" * 60)

print(f"A : {adenine_count:,}")
print(f"T : {thymine_count:,}")
print(f"G : {guanine_count:,}")
print(f"C : {cytosine_count:,}")

print("-" * 60)

print(f"GC Content : {gc_content:.2f} %")
print(f"AT Content : {at_content:.2f} %")

print(f"Nucléotide majoritaire : {most_frequent_nucleotide}")

print("=" * 60)