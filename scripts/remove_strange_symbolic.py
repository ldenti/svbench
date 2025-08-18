import sys

for line in sys.stdin:  # open(sys.argv[1]):
    if line.startswith("#"):
        print(line, end="")
    else:
        line = line.split("\t")
        if "<" in line[4] and line[4][1] == "<":
            continue
        print("\t".join(line), end="")
