import math
import logging


def Point(x, y, z):
	return (x, y, z)

def Vector(x, y, z):
	return (x, y, z)

def length(v):
	"Return length of a vector."
	sum = 0.0
	for c in v:
		sum += c * c
	return math.sqrt(sum)

def subtract(u, v):
	"Return difference between two vectors."
	x = u[0] - v[0]
	y = u[1] - v[1]
	z = u[2] - v[2]
	return Vector(x, y, z)

def dot(u, v):
	"Return dot product of two vectors."
	sum = 0.0
	for cu, cv in zip(u, v):
		sum += cu * cv
	return sum

def cross(u, v):
	"Return the cross product of two vectors."
	x = u[1] * v[2] - u[2] * v[1]
	y = u[2] * v[0] - u[0] * v[2]
	z = u[0] * v[1] - u[1] * v[0]
	return Vector(x, y, z)

def angle(v0, v1):
	"Return angle [0..pi] between two vectors."
	cosa = dot(v0, v1) / length(v0) / length(v1)
	return math.acos(cosa)

def dihedral(p0, p1, p2, p3):
	"Return angle [0..2*pi] formed by vertices p0-p1-p2-p3."
	v01 = subtract(p0, p1)
	v32 = subtract(p3, p2)
	v12 = subtract(p1, p2)
	v0 = cross(v12, v01)
	v3 = cross(v12, v32)
	# The cross product vectors are both normal to the axis
	# vector v12, so the angle between them is the dihedral
	# angle that we are looking for.  However, since "angle"
	# only returns values between 0 and pi, we need to make
	# sure we get the right sign relative to the rotation axis
	a = angle(v0, v3)
	if dot(cross(v0, v3), v12) > 0:
		a = -a
	return a

def _computeChainPhiPsi(c):
	"Compute the phi and psi angles of all amino acids"
	for i in range(1, len(c)):
		c[i][6] = _computePhi(c[i - 1], c[i])
	for i in range(0, len(c) - 1):
		c[i][7] = _computePsi(c[i], c[i + 1])

def _computePhi(prev, this):
	"Compute the phi angle of this amino acid."
	return dihedral(prev[5], this[3], this[4], this[5])

def _computePsi(this, next):
	"Compute the psi angle of this amino acid."
	return dihedral(this[3], this[4], this[5], next[3])

def readPDB(f, fname):
	"Read in PDB format input file and return a Protein instance."
	p = []
	c = []
	rSeq = None
	rChain = None
	rType = None
	bb = {}		# backbone map
	saved = False	# has this amino acid been added already?
	for line in f:
		if line[:3] == "TER":
			# Chain break
			if c:
				p.append(c)
				c = []
			continue
		if line[:4] != "ATOM":
			continue
		atomName = line[12:16].strip()
		if atomName not in [ "N", "CA", "C" ]:
			# Ignore non-backbone atoms
			logging.debug("skipping atom %s" % atomName)
			continue
		resType = line[17:20]
		resChain = line[21]
		resSeq = int(line[22:26])
		if rSeq != resSeq or rChain != resChain or rType != resType:
			rSeq = resSeq
			rType = resType
			rChain = resChain
			bb = {}
			saved = False
		x = float(line[30:38])
		y = float(line[38:46])
		z = float(line[46:54])
		bb[atomName] = Point(x, y, z)
		logging.debug("saved atom %s, len=%d", atomName, len(bb))
		if len(bb) == 3 and not saved:
			# We saw an amino acid since all three
			# keys (N, CA and C) are in our dictionary
			aa = [rSeq, rChain, rType, bb["N"], bb["CA"], bb["C"],
				None, None]
			c.append(aa)
			saved = True
	if c:
		p.append(c)
	return p

def main():
	file = open("2cyu.pdb", 'r')
	p = readPDB(file, '2cyu.pdb')

	for c in p:
		_computeChainPhiPsi(c)

	try:
		outputFile = open("output.txt", 'w')

		for chain in p:
			for res in chain:
				if not res[6]:
					outputFile.write("%d %s %s     --   %f\n" %(res[0], res[2], res[1], res[7]))
				elif not res[7]:
					outputFile.write("%d %s %s     %f   --\n" %(res[0], res[2], res[1], res[6]))
				else:
					outputFile.write("%d %s %s     %f   %f\n" %(res[0], res[2], res[1], res[6], res[7]))
				outputFile.write("\n")
		logging.info("output.txt")
		outputFile.close()

	except Exception as e:
		raise

main()
