import argparse, gzip, sys
print "Number of arguments: ", len(sys.argv)
print "The arguments are: " , str(sys.argv)
parser = argparse.ArgumentParser()
parser.add_argument('--Input',required=True)
parser.add_argument('--Output',required=True)
parser.add_argument('--Plan',required=True)
args = parser.parse_args()

planFile=open(args.Plan,'r')
if planFile.readline().rstrip() != "shadowStart\tstart\tend\tshadowEnd\tstartPOS\tendPOS":
    raise ValueError("plan file wrong header -- make sure it is the plan file made by SplitVCFwithShadow_makePlan.py")

#shadowStarts=[]
#shadowEnds=[]
#starts=[]
#ends=[]
startPOS=[]
endPOS=[]
for line in planFile:
    lineSplit=line.rstrip().split("\t")
    startPOS.append(int(lineSplit[4]))
    endPOS.append(int(lineSplit[5]))
#    shadowStarts.append(int(lineSplit[0]))
#    starts.append(int(lineSplit[1]))
#    ends.append(int(lineSplit[2]))
#    shadowEnds.append(int(lineSplit[3]))

planFile.close()

if not args.Input.endswith(".gz"):
    raise ValueError("Input argument must ends with .gz")

if not args.Output.endswith(".gz"):
    raise ValueError("Output argument must ends with .gz")

outFile=gzip.GzipFile(args.Output,'w')

if len(args.Input.split("*")) != 2:
    raise ValueError("Input argument must contain exact one * which will be replaced by chunk number")

inputFilePiece=args.Input.split("*")


#dataline=0
for fileIndex in range(len(startPOS)):
    inFile=gzip.GzipFile(inputFilePiece[0]+str(fileIndex+1)+inputFilePiece[1],'r')

    inFileHeader=[]
    line=inFile.readline()
    while not line.startswith("CHROM"):
        inFileHeader.append(line)
        line=inFile.readline()
    inFileHeader.append(line)

    if fileIndex == 0:
        for line in inFileHeader:
            outFile.write(line)


    for line in inFile:
        position=int(line.split("\t",2)[1])
        if position >= startPOS[fileIndex] and position <= endPOS[fileIndex]:
            outFile.write(line)

#    dataline -= starts[fileIndex]-shadowStarts[fileIndex]
#    for line in inFile:
#        if dataline == 89504 or dataline == 179008:
#            print(abc)
#        if dataline >= starts[fileIndex] and dataline <= ends[fileIndex]:
#            outFile.write(line)
#        dataline += 1
#    dataline -= shadowEnds[fileIndex]-ends[fileIndex]
#    outFile.write("***"+"\n")

    inFile.close()

outFile.close()
