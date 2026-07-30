"""
===========================================================
Genome Explorer
-----------------------------------------------------------
Module : report_utils.py

Ce module contient les fonctions utilitaires utilisées
pour générer les rapports scientifiques.

Aucune analyse biologique n'est réalisée ici.
Son rôle est uniquement de faciliter la mise en forme
du rapport.

Auteur : Genome Explorer Project
===========================================================
"""

from datetime import datetime


def ligne(longueur: int = 70, caractere: str = "=") -> str:
    """
    Retourne une ligne de séparation.

    Parameters
    ----------
    longueur : int
        Nombre de caractères.

    caractere : str
        Caractère utilisé.

    Returns
    -------
    str
    """

    return caractere * longueur


def titre(texte: str) -> str:
    """
    Génère un titre principal.

    Exemple
    --------
    ============================
    MON TITRE
    ============================
    """

    return (
        "\n"
        + ligne()
        + "\n"
        + texte.upper()
        + "\n"
        + ligne()
        + "\n"
    )


def sous_titre(texte: str) -> str:
    """
    Génère un sous-titre.
    """

    return (
        "\n"
        + texte
        + "\n"
        + "-" * len(texte)
        + "\n"
    )


def paragraphe(texte: str) -> str:
    """
    Retourne un paragraphe avec un saut de ligne.
    """

    return texte.strip() + "\n\n"


def information(cle: str, valeur) -> str:
    """
    Affiche une information sous forme :

    Organisme : Escherichia coli
    """

    return f"{cle:<25}: {valeur}\n"


def date_du_jour() -> str:
    """
    Retourne la date actuelle.

    Exemple
    --------
    28 juillet 2026
    """

    return datetime.now().strftime("%d/%m/%Y")


def heure_actuelle() -> str:
    """
    Retourne l'heure actuelle.
    """

    return datetime.now().strftime("%H:%M:%S")


def section(numero: int, titre_section: str) -> str:
    """
    Génère automatiquement une section numérotée.

    Exemple
    --------
    1. Introduction
    ---------------
    """

    texte = f"{numero}. {titre_section}"

    return (
        "\n"
        + texte
        + "\n"
        + "-" * len(texte)
        + "\n"
    )


def conclusion(message: str) -> str:
    """
    Met en valeur la conclusion du rapport.
    """

    return (
        "\n"
        + ligne()
        + "\nCONCLUSION\n"
        + ligne()
        + "\n"
        + message
        + "\n"
    )