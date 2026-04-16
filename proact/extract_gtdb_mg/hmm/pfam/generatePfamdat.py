import itertools as it
import glob


listpfam = []
for filename in glob.iglob('individual_hmms/*.hmm'):
    listpfam.append(filename[:-4].split('/')[1])


outfile = open('Pfam-A.hmm.dat.draft','w')



filename='Pfam-A.hmm.dat.origin'
dictionaryofpfam = {}

with open(filename,'r') as f:
    for key,group in it.groupby(f,lambda line: line.startswith('//')):
        if not key:
            group = list(group)
#	    print group
#	    print group[1].strip()
	    dictionaryofpfam[group[2].split("   ")[1].strip()] = group
            #print(group)

#print dictionaryofpfam.keys()

for pfam in listpfam:
    info = dictionaryofpfam.get(pfam)
    for line  in info :
	outfile.write(line)
    outfile.write('//\n')

outfile.close()
