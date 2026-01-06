#!/usr/bin/env python3

import sys, argparse

parser = argparse.ArgumentParser(description = "This script is a companion tool of PanCons.  It selects the n% most divergent regions.")
parser.add_argument("-i", "--input",                   help = "The BED file, output of the PanCons tool.",     required = True)
parser.add_argument("-t", "--threshold", type = float, help = "Keep the proportion of most divergent regions", default = 0.05)
parser.add_argument("-s", "--size",      type = int,   help = "Output intervals of size not less than this",   default = 100)
parser.add_argument("-d", "--distance",  type = int,   help = "Maximum size between consecutive intervals",    default = 10)
parser.add_argument("-c", "--col",       type = int,   help = "Use the n-th column as the conservation score", default = 5)
args = parser.parse_args()

# Python 0-based, human is 1-based
args.col -= 1

scores = []
with open(args.input) as input_file:
  for line in input_file:
    if line:
      l = line.strip().split()
      start = int(l[1])
      end   = int(l[2])
      size  = end - start
      score = float(l[args.col])
      scores.extend([score] * size)

scores.sort()
threshold = scores[int(len(scores) * args.threshold)]
print(f"Threshold: {threshold}", file = sys.stderr)

current_chr   = False
current_start = False
current_end   = False
with open(args.input) as input_file:
  for line in input_file:
    if line:
      l = line.strip().split()
      chr   = l[0]
      start = int(l[1])
      end   = int(l[2])
      size  = end - start
      score = float(l[args.col])
      if (current_chr) and (chr != current_chr):
        if current_end - current_start >= args.size:
          print(f"{current_chr}\t{current_start}\t{current_end}")
        current_chr   = False
        current_start = False
        current_end   = False
      if score <= threshold:
        if not current_chr:
          current_chr   = chr
          current_start = start
          current_end   = end
        elif start - current_end <= args.distance:
          current_end = end
        else:
          if current_end - current_start >= args.size:
            print(f"{current_chr}\t{current_start}\t{current_end}")
          current_chr   = chr
          current_start = start
          current_end   = end
        
if (current_chr) and (current_end - current_start >= args.size):
  print(f"{current_chr}\t{current_start}\t{current_end}")
