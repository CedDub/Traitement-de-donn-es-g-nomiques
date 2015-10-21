
import random

NT = {0:'A', 1:'C', 2:'T', 3:'G'}

ff = open("genome.fasta","w")
ff.write('>Genome\n')
G = ""
for i in range(3000000):
    k = random.randint(0,3)
    G = G+NT[k]
for i in range(0,len(G),100):
    ff.write(G[i:i+100]+'\n')
ff.close()

sizeREAD = 96

ff = open("reads500K.fasta","w")
for i in range (500000):
    k = random.randint(0,len(G)-96)
    R = ""
    e = 0
    com = ""
    for j in range(sizeREAD):
        l = random.randint(0,100)
        if l<4:
            n = random.randint(0,3)
            R = R+NT[n]
            if NT[n]!=G[k+j]:
                com = com+str(j)+":"+NT[n]+" "
                e=e+1
        else:
            R = R+G[k+j]
    ff.write('>read'+str(i)+' '+str(k)+' '+str(e)+ '('+com+')\n')
    ff.write(R+'\n')
ff.close()

