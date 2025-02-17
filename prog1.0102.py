from Bio import Seq
from difflib import SequenceMatcher

def parse_fasta(filename):
    sequences = {}
    with open(filename, 'r') as file:
        current_sequence_name = None
        current_sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # New sequence header
                if current_sequence_name is not None:
                    sequences[current_sequence_name] = current_sequence
                current_sequence_name = line[1:]
                current_sequence = ''
            else:  # Sequence data
                current_sequence += line
        # Add the last sequence
        if current_sequence_name is not None:
            sequences[current_sequence_name] = current_sequence
    return sequences

# Example usage:
filename = r'D:\Bioinformatics\human genome protein searcher\uniprotkb_proteome_UP000005640_2024_03_13.fasta.fasta'
protein_dataset = parse_fasta(filename)

# Codon table for translation
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Function to translate nucleotide sequence to protein sequence
def translate(sequence):
    aminoacidseries = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        if codon in codon_table:
            aminoacidseries += codon_table[codon]
        else:
            aminoacidseries += '?'
    return aminoacidseries


# Function to search for protein names based on their sequences
def search_protein(aminoacidseries):
    max_ratio = -1
    best_matches = []
    for protein_name, protein_aminoacidseries in protein_dataset.items():
        ratio = SequenceMatcher(None, aminoacidseries, protein_aminoacidseries).ratio()
        if ratio > max_ratio:
            max_ratio = ratio
            best_matches = [protein_name]
        elif ratio == max_ratio:
            best_matches.append(protein_name)

    if max_ratio >= 0.9:  # Adjust this threshold according to your requirement
        return best_matches
    else:
        return None

def main():
    while True:
        # Input nucleotide sequence
        sequence = input("Enter nucleotide sequence (type 'exit' to quit): ").upper()

        if sequence == 'EXIT':
            print("Exiting the program.")
            break

        # Translate sequence to protein
        protein_sequence = translate(sequence)

        # Search for protein names
        protein_names = search_protein(protein_sequence)

        # Check if protein exists
        if protein_names:
            print("Protein names:")
            for name in protein_names:
                print("-", name)
                print("Protein sequence:", protein_dataset[name])
            break
        else:
            print("No valid protein found for the given sequence.")
            suggestions = search_closest_sequence(sequence)
            if suggestions:
                print("Did you mean:")
                for suggestion in suggestions:
                    print("-", suggestion)

def search_closest_sequence(sequence):
    top_matches = []
    for protein_name, protein_aminoacidseries in protein_dataset.items():
        ratio = SequenceMatcher(None, sequence, protein_aminoacidseries).ratio()
        if len(top_matches) < 10:
            top_matches.append((protein_name, ratio))
            top_matches.sort(key=lambda x: x[1], reverse=True)  # Sort by ratio in descending order
        else:
            if ratio > top_matches[-1][1]:
                top_matches[-1] = (protein_name, ratio)
                top_matches.sort(key=lambda x: x[1], reverse=True)  # Sort by ratio in descending order

    top_matches = [match[0] for match in top_matches if match[1] >= 0.9]  # Filter matches with ratio >= 0.8
    return top_matches


if __name__ == "__main__":
    main()
