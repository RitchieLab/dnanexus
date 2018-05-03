import argparse, gzip
from math import floor
parser = argparse.ArgumentParser()
parser.add_argument('--Input',required=True)
parser.add_argument('--Output',required=True)
parser.add_argument('--WindowSizeBps',required=True)
parser.add_argument('--RoughChunkSize',required=True)
args = parser.parse_args()

roughChunkSize=int(args.RoughChunkSize)
windowSizeBps=int(args.WindowSizeBps)

inFile=gzip.GzipFile(args.Input,'r')
outFile=open(args.Output,'w')
outFile.write("\t".join(["shadowStart","start","end","shadowEnd","startPOS","endPOS"])+"\n")

line=inFile.readline()
while not line.startswith("#CHROM"):
    line=inFile.readline()

header=line.rstrip().split("\t")
posI = header.index("POS")

positions=[]
for line in inFile:
    positions.append(int(line.rstrip().split("\t",10)[posI]))

inFile.close()

chunks=len(positions)//roughChunkSize+1
print ("chunks="+str(chunks))

starts= [int(floor(i * len(positions)*1.0//chunks)) for i in range(chunks)]
ends=[i-1 for i in starts[1:chunks]]+[len(positions)-1]

shadowStarts=[]
currentPosition=0
for start in starts:
    while positions[currentPosition] < positions[start]-windowSizeBps:
        currentPosition +=1
    shadowStarts.append(max(0,currentPosition-1));

revShadowEnds=[]
currentPosition=len(positions)-1
for end in reversed(ends):
    while positions[currentPosition] > positions[end]+windowSizeBps:
        currentPosition -=1
    revShadowEnds.append(min(currentPosition+1,len(positions)-1))

shadowEnds=list(reversed(revShadowEnds))

for i in range(chunks):
    outFile.write("\t".join([str(item) for item in [shadowStarts[i],starts[i],ends[i],shadowEnds[i],positions[starts[i]], positions[ends[i]]]])+"\n")

outFile.close()

