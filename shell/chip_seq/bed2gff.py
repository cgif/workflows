import sys

##Input
if len(sys.argv)!=3:
	print("usage: python bed2gff.py inFile outFile")
	sys.exit()

##Reading file and processing each line
inFileName = sys.argv[1]
inFile = open(inFileName)

outFileName = sys.argv[2]
outFile = open(outFileName,'w')

lineList = []
count = 1

for line in inFile:
	line = line.strip().split('\t')
	chrom = line[0]
	start = line[1]
	end = line[2]
	peak_info = line[3]
	gffLine = chrom + '\t' + 'bed2gff' + '\t' + 'feature' + '\t' + str(start) + '\t' + str(end) + '\t'+ '0' + '\t' + '+' + '\t' + '.' + '\t' + 'ID, group'count+=1; 'peak_name', peak_info
	outFile.write(gffLine + '\n')

inFile.close()
outFile.close()  
	
