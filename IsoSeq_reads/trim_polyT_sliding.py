#!/usr/bin/env python3
#usage python3 trim_polyT_sliding_fastq.py input_polyA_trimmed.fastq output_polyT_trimmed.fastq


import sys

def trim_polyT_start_sliding_window(seq, window_size=10, max_mismatches=1):
    seq = seq.upper()
    n = len(seq)
    
    if n < window_size:
        mismatches = sum(1 for b in seq if b != 'T')
        if mismatches <= max_mismatches:
            return ''
        else:
            return seq
    
    trim_end = 0
    
    for i in range(n - window_size + 1):
        window = seq[i:i+window_size]
        mismatches = sum(1 for b in window if b != 'T')
        
        if mismatches > max_mismatches:
            trim_end = i
            break
    else:
        trim_end = n - window_size + 1
    
    for j in range(trim_end + window_size - 1, n):
        if seq[j] != 'T':
            trim_end = j
            break
    else:
        trim_end = n

    return seq[trim_end:]


def fastq_iter(handle):
    while True:
        name = handle.readline().rstrip()
        if not name:
            break
        seq = handle.readline().rstrip()
        plus = handle.readline().rstrip()
        qual = handle.readline().rstrip()
        yield name, seq, plus, qual


def main(fastq_in, fastq_out):
    with open(fastq_in) as fin, open(fastq_out, 'w') as fout:
        for name, seq, plus, qual in fastq_iter(fin):
            trimmed_seq = trim_polyT_start_sliding_window(seq)
            trimmed_qual = qual[len(seq) - len(trimmed_seq):]  # trim qual accordingly
            
            # Skip reads trimmed to zero length
            if len(trimmed_seq) == 0:
                continue
            
            fout.write(f"{name}\n{trimmed_seq}\n{plus}\n{trimmed_qual}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.fastq output_trimmed.fastq")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
