import argparse

def filter_fasta(input_fasta, output_fasta):
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        while True:
            header = infile.readline().strip()
            sequence = infile.readline().strip()
            
            if not header or not sequence:  # Stop when the file ends
                break
            
            # Extract the count after the last '-'
            try:
                count = int(header.split("-")[-1])
                if count > 1:
                    outfile.write(f"{header}\n{sequence}\n")
            except ValueError:
                pass  # Skip malformed headers

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter FASTA sequences with count > 1.")
    parser.add_argument("input_fasta", help="Path to the input FASTA file")
    parser.add_argument("output_fasta", help="Path to the output filtered FASTA file")

    args = parser.parse_args()
    filter_fasta(args.input_fasta, args.output_fasta)
