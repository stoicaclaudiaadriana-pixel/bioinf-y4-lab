#!/usr/bin/env python3
"""
MSA Analysis Script - Clustal Omega Only
Pentru GitHub Codespaces
"""

from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import subprocess
import os
import sys

# ============================================================================
# 0. VERIFICARE ȘI INSTALARE DEPENDINȚE
# ============================================================================

def check_and_install_tools():
    """Verifică și instalează tool-urile necesare"""
    print("="*80)
    print("VERIFICARE DEPENDINȚE")
    print("="*80)
    
    # Verifică Biopython
    try:
        import Bio
        print(f"✓ Biopython instalat (versiune {Bio.__version__})")
    except ImportError:
        print("⚠ Instalez Biopython...")
        subprocess.run([sys.executable, "-m", "pip", "install", "biopython"], check=True)
        print("✓ Biopython instalat!")
    
    # Verifică Clustal Omega
    try:
        result = subprocess.run(["clustalo", "--version"], 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE,
                              text=True)
        print(f"✓ Clustal Omega instalat")
        return True
    except FileNotFoundError:
        print("⚠ Clustal Omega nu este instalat. Instalez...")
        try:
            subprocess.run(["sudo", "apt-get", "update", "-qq"], check=True)
            subprocess.run(["sudo", "apt-get", "install", "-y", "-qq", "clustalo"], check=True)
            print("✓ Clustal Omega instalat cu succes!")
            return True
        except subprocess.CalledProcessError:
            print("✗ Nu s-a putut instala Clustal Omega automat")
            print("\nRulează manual în terminal:")
            print("  sudo apt-get update")
            print("  sudo apt-get install -y clustalo")
            return False

# Rulează verificarea
if not check_and_install_tools():
    sys.exit(1)

# ============================================================================
# 1. CITIRE ȘI SELECȚIE SECVENȚE
# ============================================================================

print("\n" + "="*80)
print("TASK 3 - MSA ONLINE (CLUSTAL OMEGA)")
print("="*80)

fasta_file = "data/work/lab/rosestoica/tp53.fa"

# Verifică dacă fișierul există
if not os.path.exists(fasta_file):
    print(f"\n⚠ ATENȚIE: Fișierul {fasta_file} nu există!")
    print("\nCaută fișiere .fa în directorul curent:")
    for root, dirs, files in os.walk(".", topdown=True):
        for file in files:
            if file.endswith('.fa') or file.endswith('.fasta'):
                print(f"  Găsit: {os.path.join(root, file)}")
    sys.exit(1)

sequences = list(SeqIO.parse(fasta_file, "fasta"))

print(f"\nTotal secvențe în dataset: {len(sequences)}")

# Selectează maxim 5 secvențe
num_sequences = min(5, len(sequences))
selected_sequences = sequences[:num_sequences]

print(f"Secvențe selectate pentru MSA: {len(selected_sequences)}")
print("\nSecvențele selectate:")
for i, seq in enumerate(selected_sequences, 1):
    print(f"  {i}. {seq.id[:60]:60} (lungime: {len(seq.seq):6} bp)")

# Salvează secvențele selectate
input_file = "selected_sequences.fasta"
SeqIO.write(selected_sequences, input_file, "fasta")
print(f"\n✓ Secvențe salvate în: {input_file}")

# ============================================================================
# 2. MSA CU CLUSTAL OMEGA
# ============================================================================

print("\n" + "="*80)
print("ALINIERE MULTIPLĂ CU CLUSTAL OMEGA")
print("="*80)

output_clustal = "alignment_clustalo.aln"
output_fasta = "alignment_clustalo.fasta"

try:
    # Rulează Clustal Omega - format Clustal
    print("\nRulează Clustal Omega...")
    clustalo_cline = ClustalOmegaCommandline(
        infile=input_file,
        outfile=output_clustal,
        outfmt="clustal",
        verbose=True,
        auto=True,
        force=True
    )
    
    print(f"Comandă: {clustalo_cline}")
    stdout, stderr = clustalo_cline()
    
    # Salvează și în format FASTA
    clustalo_fasta = ClustalOmegaCommandline(
        infile=input_file,
        outfile=output_fasta,
        outfmt="fasta",
        auto=True,
        force=True
    )
    clustalo_fasta()
    
    print(f"\n✓ Aliniere completă!")
    print(f"✓ Format Clustal: {output_clustal}")
    print(f"✓ Format FASTA:   {output_fasta}")
    
    # Afișează alinierea
    print("\n" + "="*80)
    print("REZULTAT ALINIERE MULTIPLĂ (MSA)")
    print("="*80)
    with open(output_clustal, 'r') as f:
        content = f.read()
        print(content)
    
except Exception as e:
    print(f"\n✗ Eroare la rularea Clustal Omega: {e}")
    sys.exit(1)

# ============================================================================
# 3. ANALIZĂ ALINIERE
# ============================================================================

print("\n" + "="*80)
print("ANALIZA ALINIERII")
print("="*80)

try:
    # Citește alinierea
    alignment = AlignIO.read(output_clustal, "clustal")
    
    print(f"\nNumăr secvențe aliniate: {len(alignment)}")
    print(f"Lungime aliniere: {alignment.get_alignment_length()} poziții")
    
    # Calculează statistici
    print("\nStatistici per secvență:")
    print(f"{'Secvență':<50} {'Nucleotide':>10} {'Gap-uri':>10} {'% Identitate':>15}")
    print("-"*90)
    
    for record in alignment:
        seq_str = str(record.seq)
        gaps = seq_str.count('-')
        nucleotides = len(seq_str) - gaps
        identity_pct = (nucleotides / len(seq_str)) * 100
        record_id = record.id[:48] if len(record.id) > 48 else record.id
        print(f"{record_id:<50} {nucleotides:>10} {gaps:>10} {identity_pct:>14.2f}%")
    
    # ========================================================================
    # 4. ARBORE FILOGENETIC
    # ========================================================================
    
    print("\n" + "="*80)
    print("CONSTRUCȚIE ARBORE FILOGENETIC")
    print("="*80)
    
    try:
        # Calculează matrice de distanță
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        print("\nMatrice de distanță (identity):")
        print(dm)
        
        # Construiește arborele (UPGMA)
        constructor = DistanceTreeConstructor(calculator, 'upgma')
        tree = constructor.build_tree(alignment)
        
        # Salvează arborele
        tree_file = "phylogenetic_tree.nwk"
        Phylo.write(tree, tree_file, "newick")
        print(f"\n✓ Arbore salvat în: {tree_file}")
        
        # Afișează arborele ASCII
        print("\nArbore filogenetic (reprezentare ASCII):")
        print("-"*80)
        Phylo.draw_ascii(tree)
        
    except Exception as e:
        print(f"\n⚠ Eroare la construcția arborelui: {e}")

except Exception as e:
    print(f"\n✗ Eroare la analiza alinierii: {e}")

# ============================================================================
# 5. REZUMAT
# ============================================================================

print("\n" + "="*80)
print("REZUMAT FIȘIERE GENERATE")
print("="*80)

files_generated = [
    (input_file, "Secvențe selectate (input)"),
    (output_clustal, "Aliniere MSA (format Clustal)"),
    (output_fasta, "Aliniere MSA (format FASTA)"),
    ("phylogenetic_tree.nwk", "Arbore filogenetic (Newick)"),
]

print(f"\n{'Fișier':<35} {'Status':<10} {'Dimensiune':<15} {'Descriere'}")
print("-"*95)

for filename, description in files_generated:
    if os.path.exists(filename):
        size = os.path.getsize(filename)
        size_str = f"{size} bytes" if size < 1024 else f"{size/1024:.1f} KB"
        print(f"{filename:<35} {'✓ Există':<10} {size_str:<15} {description}")

print("\n" + "="*80)
print("✓ TASK 3 COMPLET!")
print("="*80)

print("\nFișierul principal pentru raport:")
print(f"  📄 {output_clustal}")
print("\nPentru a vizualiza:")
print(f"  cat {output_clustal}")
print(f"  cat phylogenetic_tree.nwk")