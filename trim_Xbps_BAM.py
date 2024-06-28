import gzip
import sys
import pysam
import numpy as np

_, input_bam, output_bam, = sys.argv
 
def make_reference(read, MD):
	result = []
	current_number = ''


    # split MD field
	for char in MD:
		if char.isdigit():
			current_number += char
		else:
			if current_number:
				result.append(int(current_number))
				current_number = ''
			result.append(char)

	# Check if there's a trailing number
	if current_number:
		result.append(int(current_number))


    # make reference read
	ref = list(read)
	counter = 0

	for value in result:
		if isinstance(value, int):
			counter += value
		else:
			ref[counter] = value
			counter += 1


	#print('REF: ', ''.join(ref), '\t')
	return ''.join(ref)




def replace_substitutions_with_N(input_bam, output_bam):
	with pysam.AlignmentFile(input_bam, "rb") as infile, pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
		for read in infile:
			try:
				md_field = read.get_tag("MD")
				if "I" in read.cigarstring or "D" in read.cigarstring:
					continue
				if not read.is_reverse or read.is_reverse:
					read_length=len(read.seq)
					final_read = list(read.query_sequence)
					end_stop=read_length-9
					ref_seq=make_reference(read.query_sequence, md_field)
					if len(ref_seq) == len(read.query_sequence): #checks for inserts/deletions
						for index, (ref_base, read_base) in enumerate(zip(ref_seq, read.query_sequence)):
							
							if index < 9 and read_base == 'T' and ref_base == 'C':
								read.query_qualities[index]=0
								final_read[index] = 'N'
								
							elif index > end_stop and read_base == 'A' and ref_base == 'G' :
								read.query_qualities[index]=0
								final_read[index] = 'N'
								
						q=read.query_qualities #save as editing query sequence will remove quality scores
						read.query_sequence = ''.join(final_read)
						read.query_qualities = q
						#print(read.query_sequence)
						#print("")
						outfile.write(read)

			except KeyError:
				continue # these are unmapped reads which are now removed with this script

# Provide the input and output BAM file paths
input_bam_file = input_bam
output_bam_file = output_bam

# Trim the BAM file
replace_substitutions_with_N(input_bam_file, output_bam_file)