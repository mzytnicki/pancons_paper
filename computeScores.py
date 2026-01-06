import sys

annotFileName = sys.argv[1]
bedFileName   = sys.argv[2]
chrName       = sys.argv[3]
colId         = int(sys.argv[4]) - 1

inData = set()

print("Reading annotation file.", file=sys.stderr)

with open(annotFileName) as annotFile:
  for line in annotFile:
    line = line.strip().split("\t")
    if line[0] == chrName:
      start = int(line[1])
      end   = int(line[2]) - 1
      for i in range(start, end+1):
        inData.add(i)

print("  done.\nReading BED file.", file=sys.stderr)
    
with open(bedFileName) as bedFile:
  for line in bedFile:
    line  = line.strip().split("\t")
    start = int(line[1])
    end   = int(line[2]) - 1
    for i in range(start, end+1):
      if i in inData:
        print(line[colId])

print("  done.", file=sys.stderr)
