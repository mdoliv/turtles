import csv

with open("../dadi/All.cc_cm.folded.downproj.fs") as f:
    reader = csv.reader(f, delimiter=" ")
    plain_sfs = list(reader)

plain_sfs = plain_sfs[1]
plain_sfs = [float(i) for i in plain_sfs]
L = sum(plain_sfs)

theta = 175.09
