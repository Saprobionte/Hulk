"""
Hulk.py
-------
Application graphique (GUI) de bioinformatique destinée au traitement de séquences nucléotidiques
au format FASTA. Développée avec customtkinter pour l'interface et BioPython pour les opérations
sur les séquences.

Fonctionnalités principales :
    - Super consensus   : génération d'un consensus à partir de séquences paired-end (forward/reverse)
    - Décompilateur     : éclatement d'un fichier multi-FASTA en fichiers individuels
    - Compilateur       : fusion de fichiers de séquences d'un dossier en un seul fichier FASTA
    - Reverse-Complement: calcul du reverse-complément de séquences
    - Alignement        : alignement de paires de séquences via MUSCLE
    - Consensus         : génération d'une séquence consensus à partir d'un alignement

Dépendances :
    - Python 3.x
    - tkinter (bibliothèque standard)
    - customtkinter
    - BioPython (Bio.Seq, Bio.SeqIO, Bio.AlignIO, Bio.Align.Applications)
    - MUSCLE (exécutable externe requis pour les fonctions d'alignement)
"""

import sys
import os
import subprocess

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

from pathlib import Path
from io import StringIO

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

import customtkinter

customtkinter.set_appearance_mode("dark")
customtkinter.set_default_color_theme("green")


class Hulk(customtkinter.CTk):
    """
    Classe principale de l'application Hulk.

    Génère le menu principal avec les boutons de navigation vers
    chacune des fonctionnalités.
    """

    def __init__(self):
        customtkinter.CTk.__init__(self)
        self.widget_menu()
        self.wm_minsize(250, 280)

    def widget_menu(self):
        """
        Construit et place les widgets du menu principal.

        Crée les boutons suivants :
            - Aide
            - Super consensus
            - Décompilateur de séquences
            - Compilateur de séquences
            - Reverse-Complement
            - Alignement
            - Consensus
        """
        # Bouton d'aide en haut à droite
        self.but_help = customtkinter.CTkButton(self, text="?", text_color="black", command=self.wind_show_help)
        self.but_help.pack(side="top", anchor="ne")
        self.but_help.configure(height=1, width=5)

        # Super consensus
        self.but_sup_consensus = customtkinter.CTkButton(self,
            text="Super consensus", text_color="black", font=("black", 17), command=self.wind_super_consensus)
        self.but_sup_consensus.pack(pady=5)
        self.but_sup_consensus.configure(height=1, width=23)

        # Séparateur
        separator = ttk.Separator(self, orient='horizontal')
        separator.pack(fill='x', pady=15)

        # Décompilateur de séquences 
        self.but_decompilateur = customtkinter.CTkButton(self,
            text="Décompilateur de séquences", text_color="black", font=("black", 14), command=self.wind_decompilateur)
        self.but_decompilateur.pack(pady=5)
        self.but_decompilateur.configure(height=1, width=23)

        # Compilateur de séquences
        self.but_compilateur = customtkinter.CTkButton(self,
            text="Compilateur de séquences", text_color="black", font=("black", 14), command=self.wind_compilateur)
        self.but_compilateur.pack(pady=5)
        self.but_compilateur.configure(height=1, width=23)

        # Reverse-Complement
        self.but_reverse = customtkinter.CTkButton(self,
            text="Reverse // Complement", text_color="black", font=("black", 14), command=self.wind_reverseComplement)
        self.but_reverse.pack(pady=5)
        self.but_reverse.configure(height=1, width=23)

        # Alignement
        self.but_alignement = customtkinter.CTkButton(self,
            text="Alignement", text_color="black", font=("black", 14), command=self.wind_alignement)
        self.but_alignement.pack(pady=5)
        self.but_alignement.configure(height=1, width=23)

        # Consensus
        self.but_consensus = customtkinter.CTkButton(self,
            text="Consensus", text_color="black", font=("black", 14), command=self.wind_consensus)
        self.but_consensus.pack(pady=5)
        self.but_consensus.configure(height=1, width=23)

    # -------------------------------------------------------------------------
    # Fenêtres secondaires (interfaces de chaque fonctionnalité)
    # -------------------------------------------------------------------------

    def wind_show_help(self):
        """
        Ouvre une fenêtre d'aide flottante (topmost) décrivant chaque fonctionnalité
        de l'application. La fenêtre est maintenue au premier plan.
        """
        help_window = customtkinter.CTkToplevel(self)
        help_window.attributes('-topmost', 'true') 
        help_window.title("Aide")
        help_window.wm_minsize(700, 350)

        help_text = "Super consensus : Créer un consensus à partir de séquences \"paired-end\".\n \
            Les séquences \"foward\" et \"reverse\" doivent se trouver dans deux dossiers séparés,\n\
            et chaque paire de séquence à la meme position dans chacun des deux dossiers. \n\n\
            Les autres options sont les fonctions utilisées dans \"Super consensus\", mais individuellement:\n\n\
            Décompilateur de séquences : Décompile toutes les séquences contenues dans un fichier fasta,\n\
            en autant de fichiers fasta uniques de séquences présent. \n\n\
            Compilateur de séquences : Compile tout les fichiers de séquences présent \n\
            dans un dossier en un seul fichier fasta.\n\n\
            Rervese//Complement : Permet de générer la séquence comlémentaire ou la séquence reverse-complement\n\
            de toutes les séquences contenues dans un fichier fasta.\n\n\
            Alignement : Réalise l'alignement de deux séquences, avec Muscle.\n\
            Plusieurs séquences peuvent se trouver dans les deux fichiers fasta.\n\
            Un alignement sera effectué pour chaque couple de séquences\n\n\
            Consensus : Crée un consensus à partir d'un alignement."

        help_label = customtkinter.CTkLabel(help_window, text=help_text)
        help_label.pack(side=tk.LEFT)

    def wind_super_consensus(self):
        """
        Ouvre la fenêtre de configuration du workflow 'Super consensus'.

        L'utilisateur doit sélectionner :
            - Le dossier contenant les séquences forward (.seq)
            - Le dossier contenant les séquences reverse (.seq)
            - L'exécutable MUSCLE (muscle.exe)
        """
        super_consensus_window = customtkinter.CTkToplevel(self)
        super_consensus_window.wm_minsize(400, 200)
        super_consensus_window.title("Super consensus")
        super_consensus_window.attributes('-topmost', 'true')

        # Variables de liaison pour les chemins saisis
        dir1 = customtkinter.StringVar()         # Dossier forward
        dir2 = customtkinter.StringVar()         # Dossier reverse
        dir_muscle = customtkinter.StringVar()   # Chemin vers l'exécutable MUSCLE

        # Boutons d'ouverture de dossier/fichier
        self.load_dir1 = customtkinter.CTkButton(
            super_consensus_window, text="Dossier Forward", text_color="black",
            command=lambda: self.open_directory(dir1))
        self.load_dir2 = customtkinter.CTkButton(
            super_consensus_window, text="Dossier Reverse", text_color="black",
            command=lambda: self.open_directory(dir2))
        self.load_dirmuscle = customtkinter.CTkButton(
            super_consensus_window, text="muscle.exe", text_color="black",
            command=lambda: self.open_file(dir_muscle))

        # Champs d'affichage des chemins sélectionnés
        self.entry_dir1 = customtkinter.CTkEntry(super_consensus_window, textvariable=dir1, width=300)
        self.entry_dir2 = customtkinter.CTkEntry(super_consensus_window, textvariable=dir2, width=300)
        self.entry_dirmuscle = customtkinter.CTkEntry(super_consensus_window, textvariable=dir_muscle, width=300)

        # Bouton de lancement du traitement
        self.but_lanch_consensus = customtkinter.CTkButton(
            super_consensus_window, text="Lancer", text_color="black",
            command=lambda: self.super_consensus(dir1.get(), dir2.get(), dir_muscle.get()))

        # Disposition des widgets dans la fenêtre
        self.load_dir1.pack(padx=10)
        self.entry_dir1.pack(padx=10, pady=10)
        self.load_dir2.pack(padx=10)
        self.entry_dir2.pack(padx=10, pady=10)
        self.load_dirmuscle.pack(padx=10)
        self.entry_dirmuscle.pack(padx=10, pady=10)
        self.but_lanch_consensus.pack(padx=10, pady=30)

    def wind_decompilateur(self):
        """
        Ouvre la fenêtre du décompilateur de séquences.

        L'utilisateur sélectionne un fichier FASTA multi-séquences. Le traitement
        consiste à créer un fichier FASTA individuel par séquence présente dans le
        fichier source, dans un sous-dossier 'fasta_files'.
        """
        decompilateur_window = customtkinter.CTkToplevel(self)
        decompilateur_window.wm_minsize(400, 200)
        decompilateur_window.title("Décompilateur de séquences")
        decompilateur_window.attributes('-topmost', 'true')

        dir1 = customtkinter.StringVar()   # Chemin du fichier FASTA source

        self.load_dir1 = customtkinter.CTkButton(
            decompilateur_window, text="Ouvrir fasta", text_color="black",
            command=lambda: self.open_file(dir1))
        self.entry_dir1 = customtkinter.CTkEntry(decompilateur_window, textvariable=dir1, width=300)

        self.but_lanch_consensus = customtkinter.CTkButton(
            decompilateur_window, text="Lancer", text_color="black",
            command=lambda: self.decompilateur(dir1.get()))

        self.load_dir1.pack()
        self.entry_dir1.pack(pady=5)
        self.but_lanch_consensus.pack(pady=7)

    def wind_compilateur(self):
        """
        Ouvre la fenêtre du compilateur de séquences.

        L'utilisateur sélectionne un dossier et un format de fichier (.fasta, .seq,
        .fas, .txt). Tous les fichiers du format choisi présents dans le dossier
        sont fusionnés en un fichier 'sequences_compilees.fasta'.
        """
        compilateur_window = customtkinter.CTkToplevel(self)
        compilateur_window.wm_minsize(400, 200)
        compilateur_window.title("Compilateur de séquences")
        compilateur_window.attributes('-topmost', 'true')

        dir1 = customtkinter.StringVar()    # Dossier contenant les fichiers à compiler
        format = customtkinter.StringVar()  # Extension des fichiers à traiter

        self.load_dir1 = customtkinter.CTkButton(
            compilateur_window, text="Ouvrir dossier", text_color="black",
            command=lambda: self.open_directory(dir1))
        self.entry_dir1 = customtkinter.CTkEntry(compilateur_window, textvariable=dir1, width=300)

        # Liste déroulante pour choisir le format d'entrée
        list_format = customtkinter.CTkComboBox(
            compilateur_window, variable=format, state='readonly',
            values=[".fasta", ".seq", ".fas", ".txt"])
        list_format.set("Format")

        # Le fichier de sortie sera créé à la racine du dossier sélectionné
        self.but_lanch = customtkinter.CTkButton(
            compilateur_window, text="Lancer", text_color="black",
            command=lambda: self.compile_fasta(
                dir1.get(),  format.get(), False))

        self.load_dir1.pack(padx=10)
        self.entry_dir1.pack(padx=10, pady=10)
        list_format.pack(padx=10, pady=15)
        self.but_lanch.pack(padx=10, pady=30)

    def wind_reverseComplement(self):
        """
        Ouvre la fenêtre de calcul du reverse-complément.

        Le résultat est écrit dans un nouveau fichier FASTA dans le même dossier
        que le fichier source.
        """
        RC_window = customtkinter.CTkToplevel(self)
        RC_window.wm_minsize(400, 200)
        RC_window.title("Reverse Complément")
        RC_window.attributes('-topmost', 'true')

        dir1 = customtkinter.StringVar()   # Chemin du fichier FASTA source

        self.load_dir1 = customtkinter.CTkButton(
            RC_window, text="Ouvrir fasta", text_color="black",
            command=lambda: self.open_file(dir1))
        self.entry_dir1 = customtkinter.CTkEntry(RC_window, textvariable=dir1, width=300)

        self.but_lanch = customtkinter.CTkButton(
            RC_window, text="Lancer", text_color="black",
            command=lambda: self.reverse_complement(dir1.get(), "Reverse Complement"))

        self.load_dir1.pack(padx=10)
        self.entry_dir1.pack(padx=10, pady=10)
        self.but_lanch.pack(padx=10, pady=30)

    def wind_alignement(self):
        """
        Ouvre la fenêtre d'alignement de séquences.

        L'utilisateur sélectionne deux fichiers FASTA et l'exécutable MUSCLE.
        Chaque paire de séquences (une par fichier, dans l'ordre) est alignée
        individuellement. Les résultats sont écrits dans un sous-dossier "Alignement".
        """
        alignement_window = customtkinter.CTkToplevel(self)
        alignement_window.wm_minsize(400, 200)
        alignement_window.title("Alignement avec Muscle")
        alignement_window.attributes('-topmost', 'true')

        dir1 = customtkinter.StringVar()         # Fichier FASTA 1 (forward)
        dir2 = customtkinter.StringVar()         # Fichier FASTA 2 (reverse)
        dir_muscle = customtkinter.StringVar()   # Exécutable MUSCLE

        self.load_dir1 = customtkinter.CTkButton(
            alignement_window, text="Ouvrir fasta 1", text_color="black",
            command=lambda: self.open_file(dir1))
        self.load_dir2 = customtkinter.CTkButton(
            alignement_window, text="Ouvrir fasta 2", text_color="black",
            command=lambda: self.open_file(dir2))
        self.load_dirmuscle = customtkinter.CTkButton(
            alignement_window, text="muscle.exe", text_color="black",
            command=lambda: self.open_file(dir_muscle))

        self.entry_dir1 = customtkinter.CTkEntry(alignement_window, textvariable=dir1, width=300)
        self.entry_dir2 = customtkinter.CTkEntry(alignement_window, textvariable=dir2, width=300)
        self.entry_dirmuscle = customtkinter.CTkEntry(alignement_window, textvariable=dir_muscle, width=300)

        self.but_lanch_consensus = customtkinter.CTkButton(
            alignement_window, text="Lancer", text_color="black",
            command=lambda: self.alignement_solo(dir1.get(), dir2.get(), dir_muscle.get()))

        self.load_dir1.pack(padx=10)
        self.entry_dir1.pack(padx=10, pady=10)
        self.load_dir2.pack(padx=10)
        self.entry_dir2.pack(padx=10, pady=10)
        self.load_dirmuscle.pack(padx=10)
        self.entry_dirmuscle.pack(padx=10, pady=10)
        self.but_lanch_consensus.pack(padx=10, pady=30)

    def wind_consensus(self):
        """
        Ouvre la fenêtre de génération de consensus à partir d'un alignement.

        L'utilisateur sélectionne un fichier d'alignement FASTA (deux séquences alignées). 
        La séquence consensus est calculée base par base, 
        selon le code IUPAC et sauvegardée dans le même dossier sous le nom "[id_seq1]_consensus.fasta".
        """
        consensus_window = customtkinter.CTkToplevel(self)
        consensus_window.wm_minsize(400, 200)
        consensus_window.title("Consensus de séquences")
        consensus_window.attributes('-topmost', 'true')

        dir1 = customtkinter.StringVar()   # Chemin du fichier d'alignement

        self.load_dir1 = customtkinter.CTkButton(
            consensus_window, text="Ouvrir alignement", text_color="black",
            command=lambda: self.open_file(dir1))
        self.entry_dir1 = customtkinter.CTkEntry(consensus_window, textvariable=dir1, width=300)

        self.but_lanch = customtkinter.CTkButton(
            consensus_window, text="Lancer", text_color="black",
            command=lambda: self.consensus_solo(dir1.get()))

        self.load_dir1.pack(padx=10)
        self.entry_dir1.pack(padx=10, pady=10)
        self.but_lanch.pack(padx=10, pady=30)

    # -------------------------------------------------------------------------
    # Utilitaires de sélection de fichiers / dossiers
    # -------------------------------------------------------------------------

    def open_directory(self, var):
        """
        Ouvre une boîte de dialogue pour sélectionner un dossier.

        Args:
            var (customtkinter.StringVar): Variable tkinter dans laquelle le
                chemin sélectionné sera stocké.
        """
        directory = filedialog.askdirectory(initialdir='/', title="select directory")
        return var.set(directory)

    def open_file(self, var):
        """
        Ouvre une boîte de dialogue pour sélectionner un fichier.

        Args:
            var (customtkinter.StringVar): Variable tkinter dans laquelle le
                chemin sélectionné sera stocké.
        """
        file = filedialog.askopenfilename(initialdir='/', title="select file")
        return var.set(file)

    # -------------------------------------------------------------------------
    # Traitement des séquences / fichiers
    # -------------------------------------------------------------------------

    def decompilateur(self, fasta_dir):
        """
        Décompose un fichier FASTA multi-séquences en fichiers FASTA individuels.

        Crée un sous-dossier 'fasta_files' dans le même répertoire que le fichier
        source (s'il n'existe pas déjà), puis écrit un fichier FASTA par séquence,
        nommé d'après l'identifiant de la séquence.

        Args:
            fasta (str): Chemin absolu vers le fichier FASTA multi-séquences.
        """
        # Création du dossier de sortie si inexistant
        output_dir = Path(fasta_dir).parent / "fasta_files"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # Itération sur chaque séquence du fichier et écriture individuelle
        for seq in SeqIO.parse(fasta_dir, "fasta"):
            SeqIO.write(seq,output_dir / (str(seq.id) + ".fasta"), "fasta")
        return None

    def compile_fasta(self, dir, format, SC:bool):
        """
        Lit et concatène le contenu de tous les fichiers d'un format donné d'un dossier.

        Args:
            dir    (str) : Chemin absolu du dossier source.
            format (str) : Extension des fichiers à lire (ex. '.seq', '.fasta').
            SC     (bool): Si True, retourne le contenu concaténé en mémoire (usage interne
                        par super_consensus). Si False, écrit le résultat dans un fichier
                        'sequences_compilees.fasta' dans le dossier source.

        Returns:
            str | None: Contenu FASTA concaténé (SC=True) ou None après écriture fichier (SC=False).
        """
        files = sorted(Path(dir).glob(f"*{format}"))
        seq_compil = "\n".join(f.read_text() for f in files)
        
        if SC :
            return seq_compil
        else :
            with open(Path(dir) / "sequences_compilees.fasta", "w+") as new_file:
                new_file.write(seq_compil)
            return None

    def reverse_complement(self, fasta, todo):
        """
        Calcule le reverse-complément de séquences FASTA.

        Deux modes d'utilisation selon la valeur de 'todo' :
            - "superRC"           : utilisé en interne par super_consensus(). Reçoit un
                                    itérable de SeqRecord, retourne un objet parseable SeqIO.
            - "Reverse Complement": lit le fichier 'fasta', écrit le reverse-complément
                                    dans '<nom_fichier>_ReverseComplement.fasta'.

        Args:
            fasta (str | iterable): Chemin vers le fichier FASTA (modes "Complement" /
                                    "Reverse Complement") ou itérable de SeqRecord
                                    (mode "superRC").
            todo  (str): Opération à effectuer ("superRC", "Reverse Complement").

        Returns:
            SeqIO generator | None: Générateur de SeqRecord (mode "superRC") ou None
                                    (modes fichier).
        """
        if todo == "superRC":
            # Mode interne : entrée = itérable de SeqRecord, sortie = générateur SeqRecord
            out_handle = StringIO()
            SeqIO.write([rec.reverse_complement(id=rec.id, description="")
                        for rec in fasta], 
                        out_handle, "fasta")
            return SeqIO.parse(StringIO(out_handle.getvalue()), "fasta")

        elif todo == "Reverse Complement":
            # Écriture du reverse-complément dans un nouveau fichier
            SeqIO.write([rec.reverse_complement(id=rec.id, description="")
                        for rec in SeqIO.parse(fasta, "fasta")],
                        fasta.split(".fasta")[0] + "_ReverseComplement.fasta", "fasta")
            return None
        
    def alignement(self, recs, muscle_dir):
        """
        Aligne des séquences FASTA via l'exécutable MUSCLE en mode subprocess.

        Lance MUSCLE comme processus enfant, lui transmet les séquences via stdin
        au format FASTA, et récupère l'alignement depuis stdout.

        Args:
            recs       (iterable): Itérable de SeqRecord à aligner.
            muscle_dir (str)     : Chemin absolu vers l'exécutable MUSCLE.

        Returns:
            MultipleSeqAlignment: Objet d'alignement BioPython (AlignIO).
        """
        muscle_cline = MuscleCommandline(muscle_dir)

        # Lancement de MUSCLE en sous-processus avec communication via pipes
        child = subprocess.Popen(
            str(muscle_cline),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=(sys.platform != "win32"), 
        )

        # Envoi des séquences à MUSCLE via stdin
        SeqIO.write(recs, child.stdin, "fasta")
        child.stdin.close()

        # Lecture et parsing de l'alignement depuis stdout
        return AlignIO.read(child.stdout, "fasta")

    def alignement_solo(self, fasta1_dir, fasta2_dir, muscle_dir):
        """
        Aligne paire par paire les séquences de deux fichiers FASTA et sauvegarde
        les résultats individuellement.

        Les séquences des deux fichiers sont appariées dans l'ordre (zip). Pour chaque
        paire, un alignement MUSCLE est réalisé et le résultat est écrit dans le
        sous-dossier 'Alignement' (créé deux niveaux au-dessus du premier fichier),
        nommé d'après l'identifiant de la première séquence.

        Args:
            fasta1_dir     (str): Chemin vers le premier fichier FASTA.
            fasta2_dir     (str): Chemin vers le second fichier FASTA.
            muscle_dir (str): Chemin vers l'exécutable MUSCLE.
        """
        # Création du dossier de sortie si inexistant
        output_dir = Path(fasta1_dir).parent / "Alignement"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # Alignement paire par paire et sauvegarde
        for seq1, seq2 in zip(SeqIO.parse(fasta1_dir, "fasta"), SeqIO.parse(fasta2_dir, "fasta")):
            # Construction d'un mini-FASTA en mémoire pour cette paire
            seqRecord = (">" + str(seq1.id) + "\n" + str(seq1.seq) +
                         "\n" + ">" + str(seq2.id) + "\n" + str(seq2.seq))
            ali = self.alignement(SeqIO.parse(StringIO(seqRecord), "fasta"), muscle_dir)
            SeqIO.write(ali,
                        output_dir/ (str(seq1.id) + ".fasta"),
                        "fasta")
        return None

    def consensus(self, seq1, seq2):
        """
        Génère une séquence consensus base par base à partir de deux séquences alignées.

        Itère simultanément sur les deux séquences et détermine pour chaque position
        la base de consensus selon le code IUPAC (via la méthode iupac()).

        Args:
            seq1 (Seq): Première séquence alignée (BioPython Seq).
            seq2 (Seq): Deuxième séquence alignée (BioPython Seq).

        Returns:
            str: Séquence consensus sous forme de chaîne de caractères.
        """
        consensus = ""
        for base1, base2 in zip(seq1, seq2):
            consensus += self.iupac(base1, base2)
        return consensus

    def consensus_solo(self, ali_file):
        """
        Lit un fichier d'alignement à deux séquences et génère le fichier consensus.

        Le fichier de sortie est nommé '<id_seq1>_consensus.fasta' et écrit dans
        le même dossier que le fichier d'alignement source.

        Args:
            ali_file (str): Chemin vers le fichier d'alignement FASTA (2 séquences).
        """
        alignement = AlignIO.read(ali_file, "fasta")
        consensus = self.consensus(alignement[0].seq, alignement[1].seq)

        # Écriture du consensus au format FASTA
        output_path = Path(ali_file).parent / (str(alignement[0].id) + "_consensus.fasta")
        with open(output_path, "w+") as new_file:
            new_file.write(">" + alignement[0].id + "\n" + consensus)
        return None

    def super_consensus(self, forward, reverse, muscle_dir):
        """
        Workflow complet de génération de super-consensus à partir de séquences paired-end.

        Étapes réalisées :
            1. Création d'un dossier 'Super Consensus' dans le répertoire parent
               du dossier forward.
            2. Création d'un sous-dossier 'Consensus'.
            3. Pour chaque paire forward/reverse :
                a. Compilation des fichiers .seq de chaque dossier en mémoire.
                b. Calcul du reverse-complément des séquences reverse.
                c. Alignement de chaque paire forward/RC-reverse avec MUSCLE.
                d. Génération du consensus de l'alignement.
                e. Écriture du consensus dans 'Consensus/<id>.fasta'.
            4. Compilation de tous les consensus en un unique fichier
               'consensus_all.fasta' dans le dossier 'Super Consensus'.

        Args:
            forward    (str): Chemin du dossier contenant les séquences forward (.seq).
            reverse    (str): Chemin du dossier contenant les séquences reverse (.seq).
            muscle_dir (str): Chemin vers l'exécutable MUSCLE.
        """
        # Dossier de sortie principal
        dir = str(Path(forward).parent / "Super Consensus")
        os.mkdir(dir)
        os.chdir(dir)

        # Création du sous-dossier pour les consensus individuels
        if not os.path.exists("Consensus"):
            os.mkdir("Consensus")

        # Traitement paire par paire : forward + reverse-complément du reverse
        for seq1, seq2 in zip(
            SeqIO.parse(StringIO(self.compile_fasta(forward, ".seq", SC=True)), "fasta"),
            self.reverse_complement(
                SeqIO.parse(StringIO(self.compile_fasta(reverse, ".seq", SC=True)), "fasta"),
                "superRC"
            )
        ):
            # Construction du mini-FASTA en mémoire pour la paire
            seqRecordFR = (">" + str(seq1.id) + "\n" + str(seq1.seq) +
                           "\n" + ">" + str(seq2.id) + "\n" + str(seq2.seq))

            # Alignement MUSCLE de la paire
            ali = self.alignement(SeqIO.parse(StringIO(seqRecordFR), "fasta"), muscle_dir)

            # Génération et écriture du consensus
            consensus = self.consensus(ali[0].seq, ali[1].seq)
            seq_name = ali[0].id.split(".ab1")[0]   # Suppression de l'extension .ab1 éventuelle
            with open("Consensus/" + seq_name + ".fasta", "w+") as new_file:
                new_file.write(">" + seq_name + "\n" + consensus)

        # Compilation finale de tous les consensus en un seul fichier
        self.compile_fasta(Path(dir) / "Consensus", ".fasta", SC=False)
        (Path(dir) / "Consensus" / "sequences_compilees.fasta").rename(Path(dir) / "sequences_compilees.fasta") # Déplacement dans le dossier sup
        
        print("Consensus et compilations réalisés avec succès\nLes séquences se situent dans le répertoire \"Super Consensus\"")

    def iupac(self, base1, base2):
        """
        Retourne le code IUPAC correspondant à deux bases nucléotidiques alignées.

        Implémente la table des codes IUPAC pour les bases ambiguës :
            R = A ou G   (puRine)
            Y = C ou T   (pYrimidine)
            S = C ou G   (Strong)
            W = A ou T   (Weak)
            K = G ou T   (Keto)
            M = A ou C   (aMino)
            N = base inconnue ou gap

        Règles de priorité appliquées dans l'ordre :
            1. Les deux bases sont identiques           -> retourne la base
            2. N et "-"                                 -> retourne N
            3. N ou "-" + base                          -> retourne l'autre base
            4. Les deux bases sont différentes          -> retourne le code IUPAC
            5. L'une des bases est déjà un code IUPAC   -> retourne ce code (priorité base1)
            6. Cas non résolu                           -> retourne N

        Args:
            base1 (str): Première base (nucléotide ou code IUPAC).
            base2 (str): Deuxième base (nucléotide ou code IUPAC).

        Returns:
            str: Code IUPAC résultant (un caractère).
        """
        # Bases identiques
        if base1 == base2:
            return base1

        # N / -
        if set(base1 + base2) == set("N-"):
            return "N"
        elif base1 == "N" or base1 == "-":
            return base2
        elif base2 == "N" or base2 == "-":
            return base1

        # IUPAC 
        elif set(base1 + base2) == set("AG"):
            return "R"  
        elif set(base1 + base2) == set("CT"):
            return "Y"   
        elif set(base1 + base2) == set("CG"):
            return "S"   
        elif set(base1 + base2) == set("AT"):
            return "W"  
        elif set(base1 + base2) == set("GT"):
            return "K"  
        elif set(base1 + base2) == set("AC"):
            return "M"  

        # Cas où une base est déjà un code IUPAC
        else:
            if base1 in ["M", "m", "K", "k", "W", "w", "S", "s", "Y", "y", "R", "r"]:
                return base1
            if base2 in ["M", "m", "K", "k", "W", "w", "S", "s", "Y", "y", "R", "r"]:
                return base2
            else:
                return "N"   

if __name__ == "__main__":
    app = Hulk()
    app.title("Hulk")
    app.mainloop()