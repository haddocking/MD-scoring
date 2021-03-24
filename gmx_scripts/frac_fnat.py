import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser(description='Input the reference hbond.log contact file')
parser.add_argument('-r', '--reference', help='Reference structure ', required=True)
parser.add_argument('-o', '--output', help='Output file name ', required=True)
args = parser.parse_args()

with open(args.reference, 'r') as ref_fnat:
	lines1 = ref_fnat.readlines()
	output = open(args.output, 'a+' )
	output.write('#time 		% of fnat \n')
	fnat_number = float(len(lines1))
	print(fnat_number)
	for time in range (0, 100000, 500):
		try:
			with open('hbond_{}.log'.format(time), 'r') as snapshot:
				lines2 = snapshot.readlines()
				intersect = set(lines1).intersection(lines2)
			#print(intersect)
				contacts_time = float(len(intersect))
		except FileNotFoundError:
			contacts_time = 0
		perc_contacts = round ((float(contacts_time/fnat_number)), 2)
		output.write( '{}       {} \n'.format(time, perc_contacts))
	output.close()
