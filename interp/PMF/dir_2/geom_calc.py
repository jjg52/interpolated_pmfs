#!/usr/bin/env python
import sys
import numpy as np
import math
import re

class Cons:
	coords = []

def read_file1(file1):

	with open(file1, 'r') as f:
		lines = f.readlines()

	read_on = False
	Cons.coords = []

	for line in lines:
		if '!==' in line:
                	read_on = False	
		elif read_on == True:
			if '!' in line:
				pass
			else:
				Cons.coords.append(line.split())
		elif 'BAL' in line:
			read_on = True

#print(Cons.coords)

def get_distance(a,b):
	arr = []
	for x in (a, b):
		for i in Cons.coords:
			if x in i[1]:
				arr.append(i)
	x1 = arr[0]
	x2 = arr[1]
	dist = np.sqrt((float(x1[3])-float(x2[3]))**2 + (float(x1[4])-float(x2[4]))**2 + (float(x1[5])-float(x2[5]))**2)
	return dist

def get_angle(a,b,c):
	arr = []
	for x in (a,b,c):
                for i in Cons.coords:
                        if x in i[0]:
                                arr.append(i)
	x1 = arr[0]
	x2 = arr[1]
	x3 = arr[2]

	a = np.array([float(x1[3]),float(x1[4]),float(x1[5])])
	b = np.array([float(x2[3]),float(x2[4]),float(x2[5])])
	c = np.array([float(x3[3]),float(x3[4]),float(x3[5])])
	
	ba = a - b
	bc = c - b

	cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
	angle = np.arccos(cosine_angle)

	return np.degrees(angle)

def get_dihedral(a,b,c,d):
	arr = []
	for x in (a, b, c, d):
		for i in Cons.coords:
			if x in i[0]:
				arr.append(i)
	x1 = arr[0]
	x2 = arr[1]
	x3 = arr[2]
	x4 = arr[3]

	p1 = [float(x1[3]), float(x1[4]), float(x1[5])]
	p2 = [float(x2[3]), float(x2[4]), float(x2[5])]
	p3 = [float(x3[3]), float(x3[4]), float(x3[5])]
	p4 = [float(x4[3]), float(x4[4]), float(x4[5])]

	q1 = np.subtract(p2,p1) # b - a
	q2 = np.subtract(p3,p2) # c - b
	q3 = np.subtract(p4,p3) # d - c

	q1_x_q2 = np.cross(q1,q2)
	q2_x_q3 = np.cross(q2,q3)

	n1 = q1_x_q2/np.sqrt(np.dot(q1_x_q2,q1_x_q2))
	n2 = q2_x_q3/np.sqrt(np.dot(q2_x_q3,q2_x_q3))

	u1 = n2
	u3 = q2/(np.sqrt(np.dot(q2,q2)))
	u2 = np.cross(u3,u1)

	cos_theta = np.dot(n1,u1)
	sin_theta = np.dot(n1,u2)

	theta = -math.atan2(sin_theta,cos_theta)
	theta_deg = np.degrees(theta)
	if theta_deg <= 0:
		theta_deg += 360

	return theta_deg	

# Read zmat file and parse through geometries, calculating them for coordinate file

def read_file2(file2):
	with open(file2, 'r') as f:
		lines = f.readlines()

	geoms = []

	for line in lines:
		split = line.split()
		for x in split:
			if re.search("R[0-9]",x):
				nums = []
				for j in x:
					if j not in 'R':
						#print(j, end = '')
						nums.append(j)
				temp = get_distance(nums[0],nums[1])
				print("%8.3f" % temp, end=' ')
				geoms.append(temp)
#				print(nums)
			elif re.search("A[0-9]",x):
				nums = []
				neg = False
				for j in x:
					if j not in 'A':
						if j in '-':
							neg = True
						#print(j, end = '')
						else:
							nums.append(j)
				if neg == True:
					temp = get_angle(nums[0],nums[1],nums[2])*-1.0
					print("%8.3f" % temp, end=' ')
					geoms.append(temp)
				else:
					temp = get_angle(nums[0],nums[1],nums[2])
					print("%8.3f" % temp, end=' ')
					geoms.append(temp)

			elif re.search("D[0-9]",x):
				nums = []
				neg = False
				for j in x:
					if j not in 'D':
						if j in '-':
							neg = True
						#print(j, end = '')
						else:
							nums.append(j)
				if neg == True:
					temp = get_dihedral(nums[0],nums[1],nums[2],nums[3])*-1.0
					print("%8.3f" % temp, end=' ')
					geoms.append(temp)
				else:
					temp = get_dihedral(nums[0],nums[1],nums[2],nums[3])
					print("%8.3f" % temp, end=' ')
					geoms.append(temp)
			else:
				print(x, end=' ')
		print()

if __name__ == "__main__":
	read_file1(sys.argv[1])
	read_file2(sys.argv[2])
