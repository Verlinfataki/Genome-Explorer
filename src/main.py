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

from display import (afficher_cds, 
    afficher_premiers_genes, 
    afficher_resume, 
    afficher_sequence_gene, 
    afficher_composition_proteine, 
    afficher_statistiques, 
    afficher_validation_traduction,
    afficher_genes_sans_cds,
    afficher_rapport_proteine,
    afficher_comparaison,
    afficher_alignement,
    afficher_hits_blast

    )

from genes import (obtenir_infos_cds, 
    rechercher_cds, 
    genes_sans_cds
    )

from proteins import (aligner_proteines, comparer_traductions, meilleurs_hits, 
    traduire_gene, 
    analyser_composition_proteine,
    calculer_masse,
    calculer_pi,
    calculer_hydrophobicite,
    comparer_proteines,
    lancer_blast

    )


from statistic import calculer_statistiques_genes
from genome import lire_genome, compter_annotations
from tables.tables import TableGenerator
from analysis.protein_analysis import analyser_proteines
from plots.plots import PlotGenerator
import pandas as pd




# PROGRAMME PRINCIPAL

def main():
    """
    Point d'entrée du programme.
    """

    genome_record = lire_genome()

    compteur = compter_annotations(genome_record)

    afficher_resume(genome_record, compteur)

    afficher_premiers_genes(genome_record, limite=10)

    gene_statistics= calculer_statistiques_genes(genome_record)
    
    afficher_statistiques(gene_statistics)

    afficher_sequence_gene(genome_record, "thrL")

    # Traduire le gène et récupérer la protéine produite.
    proteine = traduire_gene(genome_record, "thrL")

    # Vérifier qu'une protéine a bien été obtenue.
    if proteine is not None:

        # Analyser la composition en acides aminés.
        composition = analyser_composition_proteine(proteine)

        # Afficher le résultat de l'analyse.
        afficher_composition_proteine(composition)


    # Recherche la CDS du gène thrL
    cds = rechercher_cds(genome_record, "thrL")

    # Récupère les informations de la CDS
    infos_cds = obtenir_infos_cds(cds)

    # Affiche les informations
    afficher_cds(infos_cds)

    # Traduction avec notre algorithme
    proteine = traduire_gene(genome_record, "thrL")

    # Traduction officielle
    traduction_genbank = infos_cds["translation"]

    # Comparaison
    est_valide = comparer_traductions(proteine, traduction_genbank)

    # Affichage
    afficher_validation_traduction(est_valide)

    # Recherche les gènes sans CDS
    liste = genes_sans_cds(genome_record)

    # Affichage
    afficher_genes_sans_cds(genome_record, liste)


    # Calcul de la masse moléculaire
    masse = calculer_masse(proteine)

    # Calcul du point isoélectrique
    pi = calculer_pi(proteine)

    # Calcul de l'hydrophobicité
    hydrophobicite = calculer_hydrophobicite(proteine)

    afficher_rapport_proteine(
        infos_cds,
        proteine,
        composition,
        masse,
        pi,
        hydrophobicite
    )

    # Recherche des deux CDS
    cds1 = rechercher_cds(genome_record, "thrL")
    cds2 = rechercher_cds(genome_record, "thrA")

    # Extraction des informations
    # Première protéine
    infos1 = obtenir_infos_cds(cds1)

    # Deuxième protéine
    infos2 = obtenir_infos_cds(cds2)

    # Comparaison des protéines
    resultat = comparer_proteines(
        infos1["translation"],
        infos2["translation"]
    )

    # Affichage
    afficher_comparaison(resultat)


    alignement = aligner_proteines(
        infos1["translation"],
        infos2["translation"]
    )

    afficher_alignement(alignement)

    blast = lancer_blast(infos1["translation"])

    hits = meilleurs_hits(blast)

    afficher_hits_blast(hits)

    table_generator = TableGenerator()

    table_generator.gene_statistics(gene_statistics)

    table_generator.genes(genome_record)

    proteins = analyser_proteines(genome_record)

    table_generator.proteins(proteins)

    plot_generator = PlotGenerator()

    genes = pd.read_csv(
        "outputs/tables/genes.csv"
    )

    plot_generator.gene_length_distribution(
        genes
    )


    proteins = pd.read_csv(
        "outputs/tables/proteins.csv"
    )

    plot_generator.protein_length_distribution(
        proteins
    )


    plot_generator.molecular_weight_distribution(
        proteins
    )

    plot_generator.isoelectric_point_distribution(
        proteins
    )

    print(proteins.columns.tolist())
    plot_generator.hydrophobicity_distribution(
        proteins
    )

    plot_generator.strand_distribution(genes)

    

    



# POINT D'ENTRÉE
if __name__ == "__main__":
    main()