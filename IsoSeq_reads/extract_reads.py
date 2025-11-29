import sys
from Bio import SeqIO

# Function to load the read names from a file
def load_read_list(read_list_file):
    with open(read_list_file, 'r') as f:
        reads = {line.strip() for line in f}
    return reads

# Function to extract matching reads from FASTA
def extract_reads(fasta_file, read_list_file, output_file):
    read_list = load_read_list(read_list_file)
    
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            read_id = record.description.split()[0]
            if read_id in read_list:
                SeqIO.write(record, out_f, "fasta")

# Main function to handle command-line arguments
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_reads.py <fasta_file> <read_list_file> <output_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    read_list_file = sys.argv[2]
    output_file = sys.argv[3]
    
    extract_reads(fasta_file, read_list_file, output_file)
