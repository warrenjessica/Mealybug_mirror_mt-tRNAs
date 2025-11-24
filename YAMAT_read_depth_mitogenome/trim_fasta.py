import argparse

def trim_fasta(input_fasta, output_fasta):
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        while True:
            header = infile.readline().strip()
            sequence = infile.readline().strip()
            
            if not header or not sequence:  # Stop when the file ends
                break
            
            trimmed_sequence = sequence[:-3] if len(sequence) > 3 else ""  # Remove last 3 bases
            
            outfile.write(f"{header}\n{trimmed_sequence}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim the last 3 nucleotides from each sequence in a FASTA file.")
    parser.add_argument("input_fasta", help="Path to the input FASTA file")
    parser.add_argument("output_fasta", help="Path to the output trimmed FASTA file")

    args = parser.parse_args()
    trim_fasta(args.input_fasta, args.output_fasta)
