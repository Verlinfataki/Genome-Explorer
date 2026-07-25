"""
============================================================
Projet      : Genome Explorer
Auteur      : Verlin Fataki

Description
-----------
Lecture d'un fichier GenBank et exploration des annotations
biologiques qu'il contient.

Le programme affiche :
    - le nombre total d'annotations ;
    - les différents types d'annotations ;
    - leur fréquence.

============================================================
"""

# ------------------------------------------------------
# IMPORTATION DES MODULES

from display import afficher_premiers_genes, afficher_resume, afficher_sequence_gene, afficher_composition_proteine, afficher_statistiques
from proteins import traduire_gene, analyser_composition_proteine
from statistic import calculer_statistiques_genes
from genome import lire_genome, compter_annotations



# PROGRAMME PRINCIPAL

def main():
    """
    Point d'entrée du programme.
    """

    genome_record = lire_genome()

    compteur = compter_annotations(genome_record)

    afficher_resume(genome_record, compteur)

    afficher_premiers_genes(genome_record, limite=10)

    statistiques = calculer_statistiques_genes(genome_record)
    
    afficher_statistiques(statistiques)

    afficher_sequence_gene(genome_record, "thrL")

    # Traduire le gène et récupérer la protéine produite.
    proteine = traduire_gene(genome_record, "thrL")

    # Vérifier qu'une protéine a bien été obtenue.
    if proteine is not None:

        # Analyser la composition en acides aminés.
        composition = analyser_composition_proteine(proteine)

        # Afficher le résultat de l'analyse.
        afficher_composition_proteine(composition)

# POINT D'ENTRÉE

if __name__ == "__main__":
    main()