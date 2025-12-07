#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu (Lab 1): Calculul matricei de distanțe pentru secvențe FASTA

Scop:
  1) Citește un fișier multi-FASTA (≥3 secvențe)
  2) Calculează distanțele Hamming (dacă aceeași lungime) sau p-distance
  3) Afișează matricea de distanțe (triunghiul superior)
  4) Identifică perechea cea mai apropiată

Instrucțiuni:
  python ex02_distance_matrix.py --input data/work/lab/rosestoica/tp53.fa
  
  Opțiuni:
    --method hamming     : Hamming distance (doar pentru secvențe de aceeași lungime)
    --method p_distance  : p-distance (proporția pozițiilor diferite)
    --truncate          : Trunchează la lungimea minimă (pentru Hamming cu lungimi diferite)
"""
import argparse
from pathlib import Path
import sys
from Bio import SeqIO
from Bio import Align
from typing import List, Tuple
import numpy as np


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculează distanța Hamming între două secvențe de aceeași lungime.
    Returnează numărul de poziții diferite.
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"Secvențele trebuie să aibă aceeași lungime! ({len(seq1)} vs {len(seq2)})")
    
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def p_distance(seq1: str, seq2: str) -> float:
    """
    Calculează p-distance: proporția pozițiilor diferite.
    Funcționează pentru secvențe de aceeași lungime.
    Returnează valoare între 0 și 1.
    """
    
    differences = hamming_distance(seq1, seq2)
    return differences / len(seq1)


def truncate_sequences(sequences: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """
    Trunchează toate secvențele la lungimea minimă.
    Returnează lista de secvențe truncate.
    """
    min_length = min(len(seq) for _, seq in sequences)
    print(f"[*] Trunchez toate secvențele la lungimea minimă: {min_length} bp", file=sys.stderr)
    
    truncated = [(seq_id, seq[:min_length]) for seq_id, seq in sequences]
    return truncated

def p_distance_aligned(seq1: str, seq2: str) -> float:
    """
    Calculează p-distance pe secvențe aliniate cu globalxx (Biopython).
    Funcționează pentru secvențe de lungimi DIFERITE.
    
    globalxx = aliniere globală fără penalități pentru gap-uri la capete.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    
    # Configurare scoring
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    
    # Fără penalități la capete 
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
  
    # Aliniază
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]  #  Doar prima aliniere
    
    aligned_seq1 = str(best_alignment[0])
    aligned_seq2 = str(best_alignment[1])
    
    matches = 0
    total = 0
    
    for c1, c2 in zip(aligned_seq1, aligned_seq2):
        if c1 != '-' and c2 != '-': 
            total += 1
            if c1 != c2:
                matches += 1
    
    if total == 0:
        return 0.0
    
    return matches / total



def calculate_distance_matrix(sequences: List[Tuple[str, str]], method: str = "p_distance", use_alignment: bool = False) -> np.ndarray:
    """
    Calculează matricea de distanțe pentru toate perechile de secvențe.
    
    Args:
        sequences: Lista de tupluri (id, secvență)
        method: "hamming" sau "p_distance"
    
    Returns:
        Matrice numpy cu distanțele
    """
    n = len(sequences)
    distance_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i + 1, n):  # Doar triunghiul superior
            seq1 = sequences[i][1].upper()
            seq2 = sequences[j][1].upper()
            
            try:
                if method == "hamming":
                    dist = hamming_distance(seq1, seq2)
                elif method == "p_distance":
                    if use_alignment:
                        dist = p_distance_aligned(seq1, seq2)
                    else:
                        dist = p_distance(seq1, seq2)
                else:
                    raise ValueError(f"Metodă necunoscută: {method}")
                
                distance_matrix[i][j] = dist
                
            except ValueError as e:
                print(f"[ERROR] {e}", file=sys.stderr)
                raise
    
    return distance_matrix


def print_distance_matrix(distance_matrix: np.ndarray, seq_ids: List[str], method: str):
    n = len(seq_ids)
    
    print(f"MATRICE DE DISTANȚE ({method.upper()})")
    if method == "hamming":
        print("hamming distance (număr poziții diferite)")
    else:
        print("p-distance (proporția pozițiilor diferite = impartirea la lungime)")
    # Afișare doar triunghiul superior)
    for i in range(n):
        for j in range(i + 1, n):
            dist = distance_matrix[i][j]
            if method == "hamming":
                print(f"{i}.) {seq_ids[i]} <-> {j}.) {seq_ids[j]}: {int(dist)}")
            else:
                print(f"{i}.) {seq_ids[i]} <-> {j}.) {seq_ids[j]}: {dist:.6f}")
    
    print()


def find_closest_pair(distance_matrix: np.ndarray, seq_ids: List[str], method: str) -> Tuple[int, int, float]:
    """
    Găsește perechea de secvențe cea mai apropiată (distanța minimă).
    
    Returns:
        (index1, index2, distanță_minimă)
    """
    n = len(seq_ids)
    min_dist = float('inf')
    min_i, min_j = 0, 1
    
    for i in range(n):
        for j in range(i + 1, n):
            if distance_matrix[i][j] < min_dist:
                min_dist = distance_matrix[i][j]
                min_i, min_j = i, j
    
    return min_i, min_j, min_dist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Fișier FASTA de intrare")
    ap.add_argument("--method", choices=["hamming", "p_distance"], default="p_distance",
                    help="Metodă de calcul: hamming sau p_distance (default: p_distance)")
    ap.add_argument("--truncate", action="store_true",
                    help="Trunchează secvențele la lungimea minimă (pentru lungimi diferite)")
    ap.add_argument("--align", action="store_true",
                    help="Folosește aliniament globalxx pentru p-distance (pentru lungimi diferite)")
    ap.add_argument("--max_seqs", type=int, default=10,
                    help="Limitează numărul de secvențe procesate (opțional)")
    
    args = ap.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"[ERROR] Fișierul nu există: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    print(f"[*] Citesc secvențele din: {input_path}", file=sys.stderr)
    records = list(SeqIO.parse(str(input_path), "fasta"))
    
    if len(records) < 3:
        print(f"[ERROR] Sunt necesare cel puțin 3 secvențe! Găsite: {len(records)}", file=sys.stderr)
        sys.exit(1)
    
    # Limitare număr secvențe 
    if args.max_seqs and len(records) > args.max_seqs:
        records = records[:args.max_seqs]
        print(f"[*] Limitez la primele {args.max_seqs} secvențe", file=sys.stderr)
    
    # Pregătire secvente
    sequences = [(rec.id, str(rec.seq)) for rec in records]
    seq_ids = [seq_id for seq_id, _ in sequences]
    
    print(f"[*] Am încărcat {len(sequences)} secvențe", file=sys.stderr)
    
    # Verificare lungimi
    lengths = [len(seq) for _, seq in sequences]
    print(f"[*] Lungimi: min={min(lengths)}, max={max(lengths)}, medie={np.mean(lengths):.1f}", file=sys.stderr)
    
    # Verificare dacă trebuie să trunchiem sau să aliniem
    if len(set(lengths)) > 1:  # Lungimi diferite
        print(f"[!] Secvențele au lungimi diferite!", file=sys.stderr)
        
        if args.method == "hamming" and not args.truncate:
            print("[ERROR] Hamming necesită secvențe de aceeași lungime!", file=sys.stderr)
            print("[INFO] Folosește --truncate pentru a trunchia la lungimea minimă", file=sys.stderr)
            sys.exit(1)
        
        if args.method == "p_distance" and not args.align and not args.truncate:
            print("[ERROR] Pentru lungimi diferite, alege una din opțiuni:", file=sys.stderr)
            print("  (a) --align      : folosește aliniament globalxx", file=sys.stderr)
            print("  (b) --truncate   : trunchează la lungimea minimă", file=sys.stderr)
            sys.exit(1)
        
        if args.truncate:
            sequences = truncate_sequences(sequences)
            print(f"[OK] Secvențe trunchiate", file=sys.stderr)
        elif args.align:
            print(f"[OK] Voi folosi aliniament globalxx (opțiunea a)", file=sys.stderr)
    
    # Calcul matrice de distanțe
    print(f"[*] Calculez distanțele ({args.method})...", file=sys.stderr)
    distance_matrix = calculate_distance_matrix(sequences, args.method, use_alignment=args.align)
    
    # Afișare matrice
    print_distance_matrix(distance_matrix, seq_ids, args.method)
    
    # Găsire pereche cea mai apropiată
    i, j, min_dist = find_closest_pair(distance_matrix, seq_ids, args.method)
    
    print(f"{'='*80}")
    print(f"PERECHEA CEA MAI APROPIATĂ:")
    print(f"{'='*80}")
    print(f"Secvența 1: {seq_ids[i]}")
    print(f"Secvența 2: {seq_ids[j]}")
    
    if args.method == "hamming":
        print(f"Distanță Hamming: {int(min_dist)} poziții diferite")
        print(f"Similaritate: {(1 - min_dist/len(sequences[0][1])) * 100:.2f}%")
    else:
        print(f"p-distance: {min_dist:.6f}")
        print(f"Similaritate: {(1 - min_dist) * 100:.2f}%")


if __name__ == "__main__":
    main()