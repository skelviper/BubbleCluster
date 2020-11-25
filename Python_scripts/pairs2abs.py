import sys

# python in out len

chrLenDict = {}
chrlist = []
with open(sys.argv[3],'r') as chrLen:
    for line in chrLen:
        linelist= line.strip().split(sep="\t")
        chrlist.append(linelist[0])
        chrLenDict[linelist[0]] = int(linelist[2])
    
with open(sys.argv[1],'r') as inputPairs, open(sys.argv[2],'w') as absPairs:
    absPairs.write(inputPairs.readline())
    
    for line in inputPairs.readlines():
        newline = line.strip().split(sep="\t")
        if (newline[1] in chrlist and newline[3] in chrlist):
            newline[2] = str(int(float(newline[2])) + chrLenDict[newline[1]])
            newline[4] = str(int(float(newline[4])) + chrLenDict[newline[3]])
            absPairs.write("\t".join(newline)+"\n")