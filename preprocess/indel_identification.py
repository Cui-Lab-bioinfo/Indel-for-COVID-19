import pysam
import argparse


def indel_identification(input, output_by_indel, output_by_seq):
    '''
    This function was designed for extracting insertions and deletions from SAM file.
    '''

    samfile = pysam.AlignmentFile(input, "r")
    
    with open(output_by_indel, 'w') as f1:
        with open(output_by_seq, 'w') as f2:

            # Loop through each read in the SAM file
            for read in samfile:

                # Check if the read is mapped to the reference
                if not read.is_unmapped:

                    # Get the CIGAR string for the read
                    cigar = read.cigarstring

                    # Initialize variables for counting insertions and deletions
                    num_insertions = 0
                    num_deletions = 0

                    # Initialize variables for storing insertion and deletion locations and sequences
                    insertion_locations = []
                    insertion_sequences = []
                    deletion_locations = []
                    deletion_sequences = []

                    # Loop through each operation in the CIGAR string
                    ref_pos = read.reference_start
                    query_pos = 0
                    count_ins = 0
                    count_del = 0

                    for op, length in read.cigartuples:

                        # If the operation is an insertion, increment the insertion count
                        if op == 1:
                            count_ins += 1

                            # Get the inserted sequence
                            insertion_seq = read.query_sequence[query_pos : query_pos + length]

                            f1.write(read.query_name + "\tIns\t" + str(ref_pos + 1) + '\t' + str(length) + '\t' + insertion_seq + '\n')

                            # Move the query position forward
                            query_pos += length

                        # If the operation is a deletion, increment the deletion count
                        elif op == 2:
                            count_del += 1

                            f1.write(read.query_name + "\tDel\t" + str(ref_pos + 1) + '\t' + str(length) + '\n')

                            # Move the reference position forward
                            ref_pos += length

                        # If the operation is a match or mismatch, move both reference and query positions forward
                        elif op == 0 or op == 7 or op == 8:
                            ref_pos += length
                            query_pos += length
                
                # Summarize the number of insertion and deletion in each sequence.
                f2.write(read.query_name + '\t' + str(count_ins) + '\t' + str(count_del) + '\n')

    # Close the SAM file
    samfile.close()

'''
input = "/Raid5/lixiang/Indel_202208/mission/20230327_indel_identification/01_sequences.ref.sam"
output_by_indel = "/Raid5/lixiang/Indel_202208/mission/20230327_indel_identification/02_indel.info.tsv"
output_by_seq = "/Raid5/lixiang/Indel_202208/mission/20230327_indel_identification/03_num_indel.seq.tsv"
'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        usage="Getting indels from SAM file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--sam", required=True, help="SAM file for getting indels")
    parser.add_argument("--output", required=True, help="The output of indels.")
    parser.add_argument("--summary", required=True, help="Output the number of indels in concrete sequence")

    args = parser.parse_args()
    
    input = args.sam
    output_by_indel = args.output
    output_by_seq = args.summary
    
    indel_identification(input, output_by_indel, output_by_seq)
