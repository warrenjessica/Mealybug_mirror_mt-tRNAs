import sys

def extract_reads(fastq_file, id_file, output_file):
    # Read the list of IDs into a set for fast lookup
    with open(id_file, 'r') as f:
        ids_to_extract = set(line.strip() for line in f)
    
    # Open input FASTQ file and output file
    with open(fastq_file, 'r') as fastq, open(output_file, 'w') as output:
        write_entry = False
        entry = []
        
        for line in fastq:
            # Identify the start of a FASTQ entry
            if line.startswith('@'):
                read_id = line.split()[0][1:]  # Extract ID without '@' and optional description
                write_entry = read_id in ids_to_extract
                if write_entry:
                    entry = [line]  # Start a new entry
                else:
                    entry = []
            elif write_entry:
                entry.append(line)
            
            # Write the full entry to the output if complete
            if write_entry and len(entry) == 4:  # FASTQ entries have 4 lines
                output.writelines(entry)
                write_entry = False

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_reads.py <input.fastq> <ids.txt> <output.fastq>")
        sys.exit(1)
    
    fastq_file = sys.argv[1]
    id_file = sys.argv[2]
    output_file = sys.argv[3]

    extract_reads(fastq_file, id_file, output_file)
