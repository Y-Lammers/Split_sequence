# Splits a file of fasta / fastq sequences on either the primer
# or sequencing tags used (needs to be provided by the user)

# import the argparse module to handle the input commands
import argparse

# get the commandline arguments for the files and settings
parser = argparse.ArgumentParser(description = 'Split the input file')

parser.add_argument('-i', metavar='Input file', type=str,
			help='Enter the input sequence file')
parser.add_argument('-t', metavar='Tag file', type=str,
			help='Enther the tag file')
parser.add_argument('-m', metavar='Max mismatches', type=int,
			help='The maximum number of mismatches and the sequence', default=0)
args = parser.parse_args()

def extract_sequence (seq_file):
	
	# The SeqIO module is imported from Biopython to ease the parsing of fasta files
	from Bio import SeqIO
	
	# create an empty list and determine the file type (fasta / fastq)
	seq_list, file_type = [], open(seq_file, 'r').readline()[0]

	# retrieve the sequences based on the filetype
	if file_type == '>':
		seq_list = [seq for seq in SeqIO.parse(seq_file, 'fasta')]
		file_type = 'fasta'
	else:
		seq_list = [seq for seq in SeqIO.parse(seq_file, 'fastq')]
		file_type = 'fastq'

	# return the list of sequences and file type:
	return [seq_list, file_type]

def import_tags (tag_file):
	
	# parse through the comma seperated value file and save
	# the tag sequences in a list
	tag_list = [line.replace('\n','').split(',') for line in open(tag_file, 'r')]

	# return the tag list
	return tag_list

def split (sequence_list, tag_list, max_mis):
	
	# split the sequence list based on the tags
	split_dic = {}

	# parse through the set of sequences
	for seq in sequence_list:
		
		status = ''

		# check if the sequence matches any of the tags
		for tag in tag_list:

			# check if the sequence is longer than the tag, if not the sequence will be discarded
			if len(seq.seq) <= len(tag_list[0][1]):
				status = 'short'
		
			else:
				pos, mis_count = 0, 0
			
				# the tag and sequence are compared
				# the match will be discared if there are to many
				# mismatches between the tag and sequence
				while pos < len(tag[1]):
					if seq.seq[pos] != tag[1][pos]:
						mis_count += 1
						if mis_count > max_mis:
							break
					pos += 1

				# if the tag and sequence are similar they
				# are stored in the dictionary under the tag
				# dictionary key, the tag is removed from the sequence
				if mis_count <= max_mis:
					seq = seq[len(tag[1]):]
					try:
						split_dic[tag[0]].append(seq)
					except:
						split_dic[tag[0]] = [seq]
					status = 'sorted'			
					break
	
		# if no match if found for a sequence (no similar
		# tags of the sequence being too short), it is stored
		# as 'unsorted'		
		if status != 'sorted':
			try:
				split_dic['unsorted'].append(seq)
			except:
				split_dic['unsorted'] = [seq]

	# return the dictionary
	return split_dic

def write_seqs (sorted_seq_dic, file_format, file_dir):
	# import the Biopython module for sequence handling
	from Bio import SeqIO
	import os.path
	
	out_direc = os.path.split(file_dir)[0]

	# for each tag write a seperate file
	for tag in sorted_seq_dic:

		out_file = open(out_direc + '/' + tag + '.' + file_format, 'w')
		SeqIO.write(sorted_seq_dic[tag], out_file, file_format)
		out_file.close()

def main ():	
	
	# get the sequence dictionary and file type
	seq_type = extract_sequence(args.i)

	# retrieve the tag list
	tags = import_tags (args.t)

	# split the sequence file based on the tags
	split_dic = split(seq_type[0], tags, args.m)
	
	# write the split results to seperate files
	write_seqs (split_dic, seq_type[1], args.i)

if __name__ == "__main__":
    main()

