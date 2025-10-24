"""
BioInformatics Midterm Problem 2
Ritika Lama 
-------------------------
This program calculates the smallest number of transversion and 
transition that a set of nucleotide, from different species descended from same ancetral
DNA, went through .

Flow of the code:
1. Reads a DNA sequence from a text file, removes whitespace, newlines and species name, and converts
   all characters to uppercase.
2. Validates that the sequence contains only valid nucleotides (A, C, G, T).
3. Calculates the number of transversion, transition, and its ratio.

Variables defined:
- `base`: The nucleotides in DNA sequence.
- `filename`: The name of the file.
- `dna_sequence`: The DNA sequence read from the file.
- `current_sequence`: The position while reading the sequence.
- `purine`: An array of A and G.
- `total_transition`: The total number of transition in a set of DNA sequence.
- `pyrimidine`: An array of C and T.
- `total_transversion`: The total number of transversion in a set of DNA sequence.
- `ratio`: transversion to transition ratio.

"""

"""
    This function validates that the sequence contains only valid DNA bases(A,T,G,C).

    Parameters:
        dna_sequence (str): DNA sequence to validate that was read from the text file.

    Returns:
        boolean: False if the sequence contains any nucleotide except A, C, G, T. True otherwise.
"""
def is_valid_dna(dna_sequence):
    for base in dna_sequence:
        if base not in "ATCG":
            return False
    return True

"""
    This function reads the sequence from the file and removes any whiteline or species name making sure that only DNA sequence
    returned.

    Parameters:
        dna_sequences (list): An array of DNA sequences read from the text file.
        current_sequence (str): The position while reading the sequence.

    Returns:
        dna_sequences : An array of DNA sequences read from the text file.
"""
def read_file(filename):
    dna_sequences = []
    current_sequence = ""

    #try and except block to handle file not found error
    try:
        file = open(filename, "r")
        for line in file:
            line = line.strip()
            if line == "":
                continue
            # removes the species name, removes white spaces, and changes all the characters(nucleotides) to uppercase
            if line[0] == ">":  
                if current_sequence != "":
                    dna_sequences.append(current_sequence)
                    current_sequence = ""
            else:
                current_sequence = current_sequence + line.upper()
        if current_sequence != "":
            dna_sequences.append(current_sequence)
        file.close()
        print("Read dna sequence: ", dna_sequences)
    except:
        print("File not found, check the file name and directory.")
        return []
    #returns an array of dna sequences
    return dna_sequences


"""
    This function computes the number of transversion, transition and its ratio from a set of DNA sequences.

    Parameters:
         purines (list): An array of A and G.
         pyrimidines (list): An array of C and T.
         total_transition: The total number of transition in a set of DNA sequence.
         total_transversion: The total number of transversion in a set of DNA sequence.
         sequence_length: The length of the DNA sequence.

    Returns:
        total_transversions: The total number of transversion in a set of DNA sequence.
        total_transitions: The total number of transition in a set of DNA sequence.
        ratio: transversion to transition ratio.
"""
def compute_transversion_ratio(dna_sequences):
    # the set of purines and pyrimidines
    purines = set(['A', 'G'])
    pyrimidines = set(['C', 'T'])
    total_transitions = 0
    total_transversions = 0
    sequence_length = len(dna_sequences[0])

    #loops through each position of the DNA sequences
    for position in range(sequence_length):
        unique_bases = set(sequence[position] for sequence in dna_sequences)
        # if there is all same nucleotide
        if len(unique_bases) == 1:
            total_transitions += 0
            total_transversions += 0
        # if there are two different nucleotides but lies in the same base, increase transition by 1 if purines, or transversion by 1 if pyrimidines
        elif len(unique_bases) == 2:
            base_one, base_two = list(unique_bases)
            if (base_one in purines and base_two in purines) or (base_one in pyrimidines and base_two in pyrimidines):
                total_transitions += 1
            else:
                total_transversions += 1
        else:  #if there are more than two types of nucleotides
            has_purine = any(base in purines for base in unique_bases)
            has_pyrimidine = any(base in pyrimidines for base in unique_bases)
            total_transitions += 1
            if has_purine and has_pyrimidine:
                total_transversions += 1

    # calculates the ratio of transversion to transition if there are any transitions else sets the ratio to 0
    if total_transitions > 0:
        ratio = round(total_transversions / total_transitions, 2)
    else: 
        ratio = 0
        print("No transitions found.")
    return total_transversions, total_transitions, ratio


"""
    This is the main program workflow.

    The steps followed in this function by calling other functions are as follows:
        1. Read DNA sequence from file (CS Midterm Problem 2.txt).
        2. Validate that the sequence is a proper DNA sequence.
        3. Compute the transversion ratio.
        4. Print results.
"""
def main():
    filename = "CS Midterm Problem 2.txt"
    dna_sequences = read_file(filename)
    # Validate all sequences
    for sequence in dna_sequences:
        if not is_valid_dna(sequence):
            print("Invalid nucleotide found in sequence:", sequence)
            return  
    # calculates the transversion ratio, transition and transversion count
    transversion_count, transition_count, ratio = compute_transversion_ratio(dna_sequences)

    #print the output for the transversion count, transition count and ratio
    print("Total Transversions (V):", transversion_count)
    print("Total Transitions (S):", transition_count)
    print("Transversion Ratio (V/S):", ratio)

if __name__ == "__main__":
    main()