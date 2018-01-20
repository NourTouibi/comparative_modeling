from math import sqrt

def distanceTuple(a,b,c,d):
	if a=="-" or b=="-" or c=="-" or d=="-":
		return "#NA"
	return (distance(b,d),distance(a,c),distance(b,c),distance(a,d))

def distance(atom_a, atom_b):
	a = get_coord(atom_a)
	b = get_coord(atom_b)
	xa = a[0]
	ya = a[1]
	za = a[2]
	xb = b[0]
	yb = b[1]
	zb = b[2]
	return sqrt(((xb -xa)**2)+((yb -ya)**2)+((zb-za)**2))

def get_coord(el):
	return el.get_coord()

def getPrev(seq, index, gap):
	nindex=index-gap
	while nindex<0 or seq[nindex] =='-':
		nindex+=1
	return nindex

def getNext(seq, index, gap):
	nindex=index+gap
	while nindex>len(seq)-1 or seq[nindex] =='-':
		nindex-=1
	return nindex

#Distance = sqrt(((xb -xa)**2)+((yb -ya)**2)+((zb-za)**2))

def get_distance(atoms, start, end):
	a = atoms[getPrev(atoms,start,6)]
	b = atoms[getPrev(atoms,start,5)]
	c = atoms[getNext(atoms,end,5)]
	d = atoms[getNext(atoms,end,6)]
	return distanceTuple(a,b,c,d)

def get_distance_tuples(atoms_with_gaps):
	index = 0
	dists={}
	while index<len(atoms_with_gaps):
		aa = atoms_with_gaps[index]
		if aa == '-':
			a="-"
			b="-"
			c="-"
			d="-"
			if index-1>0 and atoms_with_gaps[index-1] != '-':
				a = atoms_with_gaps[getPrev(atoms_with_gaps,index,6)]
				b = atoms_with_gaps[getPrev(atoms_with_gaps,index,5)]
			if index+1<len(atoms_with_gaps) and atoms_with_gaps[index+1] != '-':
				c = atoms_with_gaps[getNext(atoms_with_gaps,index,5)]
				d = atoms_with_gaps[getNext(atoms_with_gaps,index,6)]
			dists[index] = distanceTuple(a,b,c,d)
		index+=1
	return dists
