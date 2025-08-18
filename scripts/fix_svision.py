import sys

for line in sys.stdin:  # open(sys.argv[1]):
    if line.startswith("#"):
        print(line, end="")
    else:
        line = line.split("\t")
        if "SVTYPE=DEL;" in line[7]:
            line[4] = "<DEL>"
        elif "SVTYPE=INS;" in line[7]:
            line[4] = "<INS>"
        else:
            continue
        print("\t".join(line), end="")
