import sys

# Script to generate a range of num zmatrices from two minima zmatrices

with open(sys.argv[1], 'r') as f:
	lines1 = f.readlines()

with open(sys.argv[2], 'r') as f:
        lines2 = f.readlines()


arrs = ['one', 'two']
zmats = [lines1, lines2]
for i in range(len(zmats)):
	arrs[i] = []
	for line in zmats[i]:
		split = line.split()
		for x in split:
			if '.' in x:
				arrs[i].append(x)

file_start = 'zma_'
# number of interpolated structures
num = 13
for j in range(num+3):
	temp = []
	for i in range(len(arrs[1])):
		temp.append((float(arrs[1][i])-float(arrs[0][i]))/(num)*(j-1) + float(arrs[0][i]))
	print(j-1)
	q = 0
	with open(file_start + str(j-1), 'w+') as f:
		for i in lines1:
			split = i.split()
			for x in split:
				if '.' in x:
					f.write(str(temp[q]) + '  ')
					q += 1
				else:
					f.write(x + '  ')
			f.write('\n')
