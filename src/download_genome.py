"""
Mission 2 : Télécharger un génome depuis NCBI

Objectif :
Télécharger automatiquement le génome d'Escherichia coli
au format FASTA et GenBank.

Auteur : Verlin Fataki
Projet : Genome Explorer
"""

# ------------------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------------------

# Entrez permet de communiquer avec les serveurs du NCBI.
# SeqIO permettra plus tard de lire les fichiers biologiques téléchargés.
from Bio import Entrez, SeqIO

# pathlib facilite la manipulation des dossiers et des fichiers,
# quel que soit le système d'exploitation (Windows, Linux ou macOS).
from pathlib import Path

# Toujours renseigner une adresse e-mail valide.
# NCBI l'utilise uniquement pour vous contacter en cas de problème lié à l'utilisation de l'API.
Entrez.email = "osoverlin742@gmail.com"


# ------------------------------------------------------------------
# IDENTIFIANT DE RÉFÉRENCE DU GÉNOME
# ------------------------------------------------------------------

# Accession du chromosome de référence d'Escherichia coli K-12 MG1655.
GENOME_ACCESSION = "U00096.3"

# ------------------------------------------------------
# DOSSIER DES DONNÉES
# ------------------------------------------------------

# Création d'un objet Path représentant le dossier "data/raw".
DATA_DIRECTORY = Path("data") / "raw"

# Crée automatiquement le dossier s'il n'existe pas encore.
DATA_DIRECTORY.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# TÉLÉCHARGEMENT DU GÉNOME AU FORMAT FASTA
# ------------------------------------------------------------------
print("=" * 60)
print("TÉLÉCHARGEMENT DU FICHIER FASTA")
print("=" * 60)

# Demander au NCBI le fichier FASTA correspondant
# à l'identifiant du génome.
fasta_handle = Entrez.efetch(
    db = "nucleotide",
    id = GENOME_ACCESSION,
    rettype = "fasta",
    retmode = "text"
)

# Nom du fichier qui sera enregistré sur l'ordinateur.
fasta_file = DATA_DIRECTORY / "ecoli_genome.fasta"

# Ouvrir un nouveau fichier en mode écriture.
with open(fasta_file, "w") as file:
    file.write(fasta_handle.read())

fasta_handle.close()

print(f"Fichier FASTA enregistré : {fasta_file}")   

# ------------------------------------------------------
# TÉLÉCHARGEMENT DU GENBANK
# ------------------------------------------------------

print("\n" + "=" * 60)
print("TÉLÉCHARGEMENT DU FICHIER GENBANK")
print("=" * 60)

genbank_handle = Entrez.efetch(
    db="nucleotide",
    id=GENOME_ACCESSION,
    rettype="gb",
    retmode="text"
)

genbank_file = DATA_DIRECTORY / "ecoli_genome.gb"

with open(genbank_file, "w") as file:
    file.write(genbank_handle.read())

genbank_handle.close()

print(f"Fichier GenBank enregistré : {genbank_file}")

print("\n" + "=" * 60)
print("Téléchargement terminé avec succès.")
print("=" * 60)