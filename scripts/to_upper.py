import sys

for line in sys.stdin:  # open(sys.argv[1]):
    if line.startswith("#"):
        print(line, end="")
    else:
        line = line.split("\t")
        line[3] = line[3].upper()
        line[4] = line[4].upper()
        print("\t".join(line), end="")
