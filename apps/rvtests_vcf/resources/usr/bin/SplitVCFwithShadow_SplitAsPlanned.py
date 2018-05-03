import argparse, gzip
parser = argparse.ArgumentParser()
parser.add_argument('--Input',required=True)
parser.add_argument('--Output',required=True)
parser.add_argument('--Plan',required=True)
args = parser.parse_args()

planFile=open(args.Plan,'r')
if planFile.readline().rstrip() != "shadowStart\tstart\tend\tshadowEnd\tstartPOS\tendPOS":
    raise ValueError("plan file wrong header -- make sure it is the plan file made by SplitVCFwithShadow_makePlan.py")

shadowStarts=[]
shadowEnds=[]
for line in planFile:
    lineSplit=line.rstrip().split("\t")
    shadowStarts.append(int(lineSplit[0]))
    shadowEnds.append(int(lineSplit[3]))

planFile.close()

if not args.Input.endswith(".gz"):
    raise ValueError("Input argument must ends with .gz")

inFile = gzip.GzipFile(args.Input,'r')
inFileHeader=[]
line=inFile.readline()
while not line.startswith("#CHROM"):
    inFileHeader.append(line)
    line=inFile.readline()

inFileHeader.append(line)

if not args.Output.endswith(".gz"):
    raise ValueError("Output argument must ends with .gz")

if len(args.Output.split("*")) != 2:
    raise ValueError("Output argument must contain exact one * which will be replaced by chunk number") 

outFilePiece=args.Output.split("*")

outFiles=[]
for fileIndex in range(len(shadowStarts)):
    outFiles.append(gzip.GzipFile(outFilePiece[0]+str(fileIndex+1)+outFilePiece[1],'w'))
    for line in inFileHeader:
        outFiles[fileIndex].write(line)

dataline=0
for line in inFile:
    for fileIndex in range(len(shadowStarts)):
        if dataline >= shadowStarts[fileIndex] and dataline<= shadowEnds[fileIndex]:
            outFiles[fileIndex].write(line)
    dataline += 1

inFile.close()
for fileIndex in range(len(shadowStarts)):
    outFiles[fileIndex].close()
