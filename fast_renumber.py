# fast_renumber.py this script reads and renumber the compiled runs results from BPP
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    header = next(infile)
    outfile.write(header)
    for i, line in enumerate(infile):
        parts = line.strip().split('\t', 1)
        if len(parts) == 2:
            outfile.write(f"{2*(i+1)}\t{parts[1]}\n")

#Usage python3 fast_renumber.py input.txt output.txt
