#!/usr/bin/env python

import sys

class periodic:
	table = {'H' : '1',
		'C' : '6',
		'N' : '7',
		'O' : '8',
		'Br' : '35'} 

with open(sys.argv[1]) as f:

	dash = []

#work out number of structures from number of dashed lines
	lines = list(f.readlines())
	for i in lines:
		split = i.split()
		if len(split) == 4 :
			if split[0] in periodic.table.keys() :
				dash.append(split)			
	
#write the cooordinate files

filename = sys.argv[1].strip('.gjf') + '.crd'

with open(filename, 'w+') as f:
	f.write('!===========================================================================\n')
	f.write('%s 1 1 ! # of atoms, residues and subsystems.\n' % len(dash))
	f.write('!===========================================================================\n')
	f.write('Subsystem	1	SOLUTE\n')
	f.write('	1 ! # of residues.\n')
	f.write('!===========================================================================\n')
	f.write('Residue	1	BAL\n')
	f.write('	%s ! # of atoms.\n' % len(dash))
	num = 1
	for j in dash:
#		split_j = str(j).split()	
		f.write('%s	%s%s	%s	%s	%s	%s\n' % (num, j[0], num, periodic.table[j[0]], j[1], j[2], j[3]))
		num += 1
	f.write('\n!===========================================================================')

with open('num_files', 'w+') as f:
	f.write('%d' % (len(dash) + 1) )

with open('num_atoms', 'w+') as f:
	f.write('%d' % len(dash))
