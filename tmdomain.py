"""
Bioinformatics Midterm Problem 1
Ritika Lama
-------------------------
This program predicts transmembrane domains in protein sequences using
the Kyte-Doolittle hydrophobicity scale. Transmembrane regions are sequence of 
18-21 amino acids that whose hydrophobicity is greater than defined threshold.

Flow of the code:
1. Reads amino acid sequences from a text file ("CS Midterm Problem 1.txt").
2. Validates that each sequence contains only valid amino acids.
3. Calculates hydrophobicity for each residue using the K-D scale.
4. Computes average hydrophobicity of a region and moves on to next until the end.
5. Identifies potential TM regions where the average hydrophobicity > 0.
6. Highlights predicted TM domains in green.
7. Prints positions and average hydrophobicity scores for each TM region.

Variables defined:
- KD: Dictionary of amino acids and their KD hydrophobicity values.
- sequences: List of amino acid sequences read from the file.
- transmembrane_region: List of tuples containing TM start, end, and average hydrophobicity.
"""

import sys
"""
    This list defines the KD hydrophobicity scale for amino acids.
    Positive values indicate hydrophobic amino acids while
    negative values indicate hydrophilic amino acids.
"""
KD = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

"""
    Function: read_sequences
    Purpose:
        Reads amino acid sequences from a text file and formats them in uppercase.
    Parameters:
        filename (str): Name of the file containing AA sequences.
    Returns:
        list: List of AA sequences (uppercase, whitespace removed).
        Prints an error message if the file is not found.
"""
def read_sequences(filename):
    try:
        with open(filename, "r") as f:
            return [line.strip().upper() for line in f if line.strip()]
    except FileNotFoundError:
        print("File not found, check the file name and directory.")
        sys.exit()

"""
    Function: validate_sequence
    Purpose:
         Validates that a given protein sequence contains only valid set of amino acids.
    Parameters:
         sequence (str): Protein sequence to validate.
         valid_amino_acids : Valid amino acid symbols from KD scale.
    Returns:
        Prints an error and halts the program if an invalid residue is found.
"""
def validate_sequence(sequence, valid_amino_acids):
    for amino_acid in sequence:
        if amino_acid not in valid_amino_acids:
            print(f"Invalid amino acid '{amino_acid}' found in sequence: {sequence}")
            sys.exit("The file contains invalid AA sequence. ")

"""
     Function: average_hydrophobicity
     Purpose:
         Calculates the mean hydrophobicity of an amino acid region.
     Parameters:
         sequence_window (str): Substring of amino acids.
         scale (dictionary): Hydrophobicity scale.
     Returns:
         float: Average hydrophobicity value for the window.
"""
def average_hydrophobicity(sequence_window, scale):
    total_hydrophobicity = 0
    # total hydrophobicity values for each amino acid in the window
    for amino_acid in sequence_window:
        total_hydrophobicity += scale.get(amino_acid, 0)
    return total_hydrophobicity / len(sequence_window)

"""
     Function: find_transmembrane_regions
     Purpose:
         Identifies regions in a protein sequence with hydrophobicity above threshold.
     Parameters:
         sequence (str): Amino acid sequence.
         scale (dictionary): Hydrophobicity scale.
         window_size (int): Number of amino acids per analysis window.
         threshold (float): Hydrophobicity > threshold(0) indicates TM region.
     Returns:
        list: Tuples of (start, end, average hydrophobicity) for TM regions.
"""
def find_transmembrane_regions(sequence, scale, window_size=19, threshold=0):
    # after calculating the hydrophobicity for every window of 19 amino acids, it returns the transmembrane region where the average hydrophobicity is greater than 0
    transmembrane_region = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        avg_hydrophobicity = average_hydrophobicity(window, scale)
        if avg_hydrophobicity > threshold:
            transmembrane_region.append((i + 1, i + window_size, avg_hydrophobicity))
    return transmembrane_region

"""
     Function: highlighted_transmembrane
     Purpose:
         Highlights predicted TM regions in green for visual output.
     Parameters:
         sequence (str): Full amino acid sequence.
         transmembrane_region (list): TM region positions.
     Returns:
         str: Sequence string with TM region highlighted in green.
"""
def highlighted_transmembrane(sequence, transmembrane_region):
    highlighted = ""
    # returns the amino acid sequence with the transmembrane region highlighted in green
    for i, amino_acid in enumerate(sequence, start=1):
        if any(start <= i <= end for start, end, _ in transmembrane_region):
            highlighted += "\033[92m" + amino_acid + "\033[0m"
        else:
            highlighted += amino_acid
    return highlighted

"""
     Function: analyze_sequences
     Purpose:
        Analyzes file, validation, hydrophobicity analysis,
         and visualization for all sequences in the file.
     Parameters:
         filename (str): Path to the file containing amino acid sequences.
"""
def analyze_sequences(filename):
    sequences = read_sequences(filename)

    # adds a counter to sequences and returns an enumerate object.
    for index, sequence in enumerate(sequences, start=1):
        #validate sequence, checks if the sequence contains only valid amino acid same as the keys for the KD scale
        validate_sequence(sequence, KD.keys())
        print(f"\nSequence {index}: {sequence}")
        transmembrane_region = find_transmembrane_regions(sequence, KD)

        if transmembrane_region:
            print("Predicted TM regions (highlighted in green):")
            print(highlighted_transmembrane(sequence, transmembrane_region))
            #Prints out average hydrophobicty for each window
            print("TM region positions and hydrophobicity scores:")
            for start, end, score in transmembrane_region:
                print(f"  {start}-{end} -> Average hydrophobicity: {score:.2f}")
        else:
            print("There are no transmembrane regions detected.")

"""
     Function: main
     Purpose:
         Main function that starts TM domain analysis.
     Steps:
         1. Reads amino acid sequences from the file.
         2. Validates input sequences.
         3. Predicts TM domains using KD scale.
         4. Displays results in the terminal.
"""
def main():
    analyze_sequences("CS Midterm Problem 1.txt")

# main execution
if __name__ == "__main__":
    main()
