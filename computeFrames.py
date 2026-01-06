import sys

inputFileName = sys.argv[1]
outputPrefix  = sys.argv[2]

outputFileNames = [f"{outputPrefix}{s}.bed" for s in range(3)]
outputFiles     = [open(outputFileName, 'w') for outputFileName in outputFileNames]

with open(inputFileName) as inputFile:
  for line in inputFile:
    if line:
      if line[0] != '#':
        line = line.strip().split("\t")
        if line[2] == "CDS":
          ref    = line[0]
          start  = int(line[3])
          end    = int(line[4])
          frame  = int(line[7])
          strand = line[6]
          if strand == '+':
            pos = range(start, end + 1)
          else:
            pos = range(end, start - 1, -1)
          for i in pos:
            outputFiles[frame].write(f"{ref}\t{i-1}\t{i}\n")
            frame = (frame + 1) % 3
