#!/usr/bin/env python
'''calculate mutual information of two lists of numbers or symbols'''

import numpy as np

from collections import Counter

__author__ = "Liguo Wang"
__copyright__ = ""
__credits__ = []
__license__ = "GPLv2"
__version__ = "1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "Wang.Liguo@mayo.edu"
__status__ = "Development" #Prototype or Production

def Mutual_information1(x,y):
	'''
	x and y are lists of symbols (like 'A','C','G','T').
	Calculation mutual information based on: MI = H(x) + H(y) - H(x,y)
	Log2 based, unit is bit
	'''
	x = [str(i) for i in x]
	y = [str(i) for i in y]
	if len(x) != len(y):
		return 0
	xy = [''.join(i) for i in zip(x,y)]
	
	x_freq = np.array(list(Counter(x).values()))	#a.items(): [('A', 3), ('C', 2), ('T', 4), ('G', 1)]
	y_freq = np.array(list(Counter(y).values()))
	xy_freq = np.array(list(Counter(xy).values()))
	
	x_freq = x_freq*1.0/sum(x_freq)
	y_freq = y_freq*1.0/sum(y_freq)
	xy_freq = xy_freq*1.0/sum(xy_freq)
	
	x_H = -sum([i * np.log2(i) for i in x_freq])
	y_H = -sum([i * np.log2(i) for i in y_freq])
	xy_H = -sum([i * np.log2(i) for i in xy_freq])
	
	return (x_H,y_H,xy_H, x_H+y_H-xy_H)

def Mutual_information2(x,y):
	'''
	x and y are lists of symbols (like 'A','C','G','T'). 
	Calculate mutual information based on its original definition. 
	Log2 based, unit is bit
	'''
	x = [str(i) for i in x]
	y = [str(i) for i in y]
	if len(x) != len(y):
		return 0
	xy = [''.join(i) for i in zip(x,y)]
	
	px = {}
	py = {}
	pxy = {}
	for i,j in list(Counter(x).items()):
		px[i] = j*1.0/len(x)
	for i,j in list(Counter(y).items()):
		py[i] = j*1.0/len(y)
	for i,j in list(Counter(xy).items()):
		pxy[i] = j*1.0/len(xy)
	#print px
	#print py
	#print pxy
	
	mi_sum = 0.0
	tmp = set()
	for xi, yi in zip(x,y):
		xyi = xi + yi
		if xyi in tmp: continue
		#print "%s::px:%f, py:%f, pxy:%f" % (xyi,px[xi],py[yi],pxy[xyi])
		mi_sum += (pxy[xyi] * np.log2(pxy[xyi] / (px[xi] * py[yi])))
		tmp.add(xyi)
	return mi_sum

def PMI(x,y):
	'''
	x and y are lists of symbols (like 'A','C','G','T'). 
	Calculate pointwise mutual information based on its original definition. 
	Log2 based, unit is bit
	'''
	x = [str(i) for i in x]
	y = [str(i) for i in y]
	
	if len(x) != len(y):
		return 0
	xy = [''.join(i) for i in zip(x,y)]
	#print xy
	#print set(x)
	#print set(y)
	px = {}
	py = {}
	pxy = {}
	for i,j in list(Counter(x).items()):
		px[i] = j*1.0/len(x)
	for i,j in list(Counter(y).items()):
		py[i] = j*1.0/len(y)
	for i,j in list(Counter(xy).items()):
		pxy[i] = j*1.0/len(xy)
	#print px
	#print py
	#print pxy
	
	#print set(x)
	#print set(y)

	for i in set(x):
		#if px[i] < 0.05:continue
		for j in set(y):
			#if py[j] < 0.05:continue
			tmp1 = i + j
			if i + j in pxy:
				tmp2 = np.log2(pxy[i+j]/(px[i] * py[j]))
				#tmp2 = -(np.log2(pxy[i+j]/(px[i] * py[j]))) / np.log2(pxy[i+j])
			else:
				continue
			print(tmp1,tmp2)
	
def Mutual_expected():
	'''
	x and y are lists of symbols (like 'A','C','G','T').
	Calculation mutual information based on: MI = H(x) + H(y) - H(x,y)
	Log2 based, unit is bit
	'''
	x_freq = [0.25]*4
	y_freq = [0.25]*4
	xy_freq= [0.0625]*16
	
	x_H = -sum([i * np.log2(i) for i in x_freq])
	y_H = -sum([i * np.log2(i) for i in y_freq])
	xy_H = -sum([i * np.log2(i) for i in xy_freq])
	
	return (x_H,y_H,xy_H, x_H+y_H-xy_H)

	
if __name__=='__main__':
	x=['G', 'T', 'C', 'A', 'T', 'T', 'A', 'C', 'T', 'A']
	y=['A', 'C', 'A', 'C', 'A', 'A', 'G', 'A', 'G', 'A']
	z=['A', 'C', 'A', 'T', 'A', 'A', 'T', 'A', 'G', 'T']
	
	
	X_0 = (0, 0, 1, 1, 0, 1, 1, 2, 2, 2)
	X_1 = (3, 4, 5, 5, 3, 2, 2, 6, 6, 1)
	X_2 = [7, 2, 1, 3, 2, 8, 9, 1, 2, 0]
	
	C8 = ['T', 'A', 'A', 'A', 'T', 'T', 'G', 'C', 'T', 'A', 'A', 'A', 'A', 'T', 'T', 'T', 'A', 'T', 'T', 'A', 'G', 'A', 'C', 'A', 'G', 'G', 'G', 'G', 'T', 'A', 'T', 'T', 'A', 'G', 'C', 'T', 'T', 'T', 'A', 'T', 'G', 'G', 'A', 'C', 'T', 'C', 'C', 'C', 'C', 'T', 'G', 'A', 'T', 'A', 'T', 'A', 'A', 'T', 'C', 'C', 'T', 'G', 'A', 'C', 'G', 'T', 'T', 'C', 'T', 'T', 'A', 'T', 'A', 'A', 'A', 'A', 'T', 'G', 'A', 'A', 'T', 'G', 'A', 'C', 'T', 'A', 'C', 'A', 'A', 'C', 'C', 'A', 'A', 'A', 'C', 'C', 'T', 'T', 'T', 'A', 'A', 'A', 'G', 'A', 'C', 'T', 'A', 'C', 'T', 'T', 'T', 'T', 'T', 'T', 'G', 'A', 'A', 'T', 'A', 'G', 'G', 'A', 'A', 'A', 'T', 'A', 'G', 'G', 'C', 'C', 'T', 'C', 'G', 'A', 'T', 'T', 'G', 'T', 'T', 'C', 'G', 'G', 'C', 'T', 'T', 'G', 'A', 'G', 'T', 'A', 'T', 'T', 'T', 'A', 'T', 'A', 'C', 'C', 'C', 'C', 'T', 'T', 'T', 'A', 'C', 'C', 'A', 'A', 'C', 'T', 'C', 'A', 'C', 'T', 'T', 'T', 'T', 'A', 'G', 'T', 'A', 'G', 'T', 'C', 'A', 'T', 'C', 'A', 'C', 'C', 'C', 'T', 'T']
	C9 = ['T', 'T', 'G', 'A', 'C', 'T', 'A', 'C', 'G', 'C', 'G', 'G', 'G', 'G', 'G', 'T', 'C', 'G', 'T', 'A', 'G', 'G', 'T', 'G', 'C', 'T', 'T', 'C', 'G', 'T', 'G', 'G', 'C', 'T', 'C', 'A', 'T', 'T', 'C', 'A', 'T', 'A', 'A', 'C', 'C', 'G', 'T', 'G', 'C', 'C', 'T', 'C', 'G', 'C', 'T', 'T', 'C', 'C', 'T', 'A', 'T', 'C', 'G', 'A', 'C', 'T', 'C', 'T', 'T', 'C', 'A', 'T', 'G', 'A', 'A', 'G', 'T', 'C', 'T', 'T', 'T', 'C', 'A', 'T', 'T', 'T', 'T', 'T', 'G', 'A', 'A', 'C', 'T', 'G', 'C', 'C', 'T', 'C', 'C', 'A', 'G', 'G', 'G', 'G', 'G', 'A', 'G', 'A', 'G', 'T', 'T', 'C', 'T', 'A', 'C', 'G', 'G', 'T', 'G', 'C', 'C', 'G', 'A', 'A', 'G', 'A', 'G', 'T', 'A', 'C', 'A', 'C', 'C', 'G', 'G', 'C', 'A', 'C', 'G', 'C', 'A', 'A', 'A', 'T', 'G', 'A', 'T', 'T', 'A', 'A', 'G', 'A', 'T', 'G', 'A', 'A', 'A', 'A', 'C', 'A', 'G', 'C', 'C', 'A', 'C', 'C', 'G', 'A', 'C', 'C', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'A', 'A', 'C', 'T', 'A', 'G', 'A', 'G', 'C', 'T', 'A', 'T', 'C', 'A', 'T', 'G']
	
	#AR motif pos-5 and pos-11
	C5 =  ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'G', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'G', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'T', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'T', 'C', 'C', 'C', 'C', 'T', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
	C11 = ['G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'C', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'A', 'T', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'A', 'G', 'G', 'G', 'G', 'G', 'A', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'A', 'G', 'G', 'G', 'G', 'G']
	#PMI(C5,C11)


	#https://en.wikipedia.org/wiki/Pointwise_mutual_information example
	#x=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]
	#y=[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1]
	#z=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0]
	#PMI(x,z)
	C1 = ['A','C','G','T']
	C2 = ['A','C','G','T']
	a=Mutual_information1(C1,C2)
	print(a)
	b=Mutual_information2(C1,C2)
	print(b)
