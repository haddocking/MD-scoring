import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser(description='Calculates energy differences from the first snapshot')
parser.add_argument('-e', '--efile', help='Energy file', required=True)
parser.add_argument('-o', '--output', help='Output file name ', required=True)
args = parser.parse_args()


with open(args.efile , 'r') as tser:
	output = open(args.output, 'w' )
	#x=(tser.readline().split(' ')[1])
	first_line = tser.readline().split()
	time0, first_val= first_line[0],  first_line[1]
	print(time0, first_val)
	output.write('{} {} \n'.format(time0, (float(first_val)-float(first_val))))
	for line in tser:
		 #print(line[15])
		time, snapshot = line.split()
		# print(time)
		dt=(float(float(snapshot)-float(first_val)))
		output.write( '{} {} \n'.format(time, dt))
		output.close
