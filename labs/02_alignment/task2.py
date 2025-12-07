from Bio import SeqIO
from Bio import Align

file_path = "data/work/lab/rosestoica/tp53.fa"
sequences = list(SeqIO.parse(file_path, "fasta"))

if len(sequences) < 2:
    print("Nu sunt suficiente secvențe!")
else:
    seq1 = sequences[0][5:500]
    seq2 = sequences[1][0:500]
    
    print(f"Secvențe: {seq1.id} vs {seq2.id}")
    print("-" * 40)

    aligner = Align.PairwiseAligner()
    
    # echivalent cu globalxx/localxx: match=1, mismatch=0, gap=0
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = 0.0
    aligner.extend_gap_score = 0.0

    # ALINIERE GLOBALĂ
    aligner.mode = 'global'
    alignments_global = aligner.align(seq1.seq, seq2.seq)
    best_global = alignments_global[0] # Luăm primul rezultat

    print("\n### 1. ALINIERE GLOBALĂ (Fragment) ###")
    print(f"Scor: {best_global.score}")
    print(best_global) 

    # ALINIERE LOCALĂ
    aligner.mode = 'local'
    alignments_local = aligner.align(seq1.seq, seq2.seq)
    best_local = alignments_local[0]

    print("\n### 2. ALINIERE LOCALĂ (Fragment) ###")
    print(f"Scor: {best_local.score}")
    print(best_local)
