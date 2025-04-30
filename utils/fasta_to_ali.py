#!/usr/bin/env python3
import argparse


def read_fasta(fp):
    name, seq = None, []
    # Read the opened file line by line
    for line in fp:
        line = line.rstrip()
        # If the line starts with a >, it is the name of the sequence
        if line.startswith(">"):
            if name:
                # If the name is not None, it means that we have already read a sequence
                yield name, ''.join(seq)
            name, seq = line[1:], []
        # Else, it is the sequence
        else:
            seq.append(line)
    if name:
        # Once we have read the last sequence, we yield it
        yield name, ''.join(seq)


def fasta_parser(fasta_file):
    names, seqs = [], []
    # Read the fasta file
    with open(fasta_file, 'r') as fp:
        # For each name and sequence in the fasta file, append them to the list
        for name, seq in read_fasta(fp):
            names.append(name)
            seqs.append(seq)
    # return the list of names and sequences
    return names, seqs


def main(args):
    # Read the input file in fasta format
    # Create a list with the names and sequences
    names, seqs = fasta_parser(args.input)
    # Assert that the number of sequences is equal to the number of names
    assert len(names) == len(seqs)

    # Assert that all sequences have the same length
    set_seq_len = set([len(seq) for seq in seqs])
    assert len(set_seq_len) == 1, "All sequences must have the same length"
    seq_len = set_seq_len.pop()

    # Open the output file
    ali_file = open(args.output, 'w')
    # Write the header (the number of sequences and the sequence length)
    ali_file.write(f"{len(seqs)} {seq_len}\n")
    # Write the names and sequences
    ali_file.write("\n".join([f"{name} {seq}" for name, seq in zip(names, seqs)]))
    ali_file.close()


if __name__ == '__main__':
    # Parse arguments
    descr = "A script to convert a fasta file to a phylip file.\n"
    descr += "The script will read the fasta file and write it in a new file in phylip format.\n\n"
    descr += "The fasta format is alternating the name of the sequence (starting with >) and the sequence itself in the next line:\n"
    descr += ">name1\nATCTC\n>name2\nATGTA\n>name3\nATGTC\n\n"
    descr += "For the phylip, the first line is the number of sequences and the length of the sequences, with then name and sequence in the same line:\n"
    descr += "3 5\nname1 ATCTC\nname2 ATGTA\nname3 ATGTC\n"

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=descr)
    parser.add_argument('-i', '--input', required=True, type=str, dest="input",
                        help="Input alignment file in fasta format")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output",
                        help="Output alignment file in phylip format")

    # Open alignment file in fasta format and write it in a new file in phylip format
    main(parser.parse_args())
