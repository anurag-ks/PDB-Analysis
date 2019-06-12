import math
import matplotlib.pyplot as plt
import Bio.PDB


parser = Bio.PDB.PDBParser()
builder = Bio.PDB.PPBuilder()
structure = parser.get_structure("2CYU", "2cyu.pdb")


def degrees(radians):
    if radians:
        angle = radians*180 / math.pi
        while angle > 180:
            angle = 360 - angle
        while angle < -180:
            angle = 360 + angle
        return angle
    return None

phi_psi_degrees = []

for model in structure:
    temp = []
    for chain in model :
        polypeptides = builder.build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides) :
            print("Model %s Chain %s" % (str(model.id), str(chain.id)),)
            print("(part %i of %i)" % (poly_index+1, len(polypeptides)),)
            print("length %i" % (len(poly)),)
            print("from %s%i" % (poly[0].resname, poly[0].id[1]),)
            print("to %s%i" % (poly[-1].resname, poly[-1].id[1]))
            phi_psi = poly.get_phi_psi_list()
            for res_index, residue in enumerate(poly) :
                res_name = "%s%i" % (residue.resname, residue.id[1])
                phi, psi = phi_psi[res_index]
                temp.append((degrees(phi), degrees(psi)))
                print(res_name,)
                print((degrees(phi), degrees(psi)))
    phi_psi_degrees.append(temp)
    temp = []

models_count = len(phi_psi_degrees)
fig = plt.figure(figsize=(50, 50))
for i in range(int(models_count)):
    plt.subplot(5, 4, i+1)
    plt.title("Model %s" % str(i))
    for i in phi_psi_degrees[i]:
        if(i[0] and i[1]):
            plt.scatter(i[0], i[1])
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\psi$')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid()
fig.subplots_adjust(hspace=1.5, wspace=0.5, top=0.95)
fig.savefig("full.png", bbox_inches='tight')
plt.close(fig)

for i in range(models_count):
    fig = plt.figure(figsize=(10, 10))
    plt.title("Model %s" %(str(i)))
    for j in phi_psi_degrees[i]:
        if(j[0] and j[1]):
            plt.scatter(j[0], j[1])
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\psi$')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid()
    plt.savefig("model%s.png"%str(i), bbox_inches='tight')
    plt.close(fig)
    print(str(i),)
