#!/usr/bin/env python3
import re
import csv

# 0) Set your filenames here:
infile   = 'matchtable_id.txt'           
remfile  = 'paralogs.txt'  
outfile  = 'without-paralogs.txt'

# 1) Define your FOâ€‘headers (87 columns):
header_line = """FO0090	FO0001	FO0028	FO0046	FO0049	FO0061	FO0066	FO0067	FO0069	FO0070
FO0079	FO0080	FO0087	FO0088	FO0089	FO0000	FO0091	FO0092	FO0094	FO0096
FO0097	FO0098	FO0099	FO0100	FO0104	FO0105	FO0107	FO0111	FO0112	FO0113
FO0118	FO0119	FO0124	FO0128	FO0129	FO0130	FO0132	FO0136	FO0138	FO0139
FO0145	FO0146	FO0149	FO0150	FO0151	FO0154	FO0156	FO0157	FO0168	FO0169
FO0171	FO0172	FO0173	FO0174	FO0176	FO0177	FO0178	FO0179	FO0180	FO0181
FO0182	FO0183	FO0184	FO0185	FO0186	FO0187	FO0188	FO0189	FO0191	FO0192
FO0203	FO0210	FO0219	FO0221	FO0224	FO0227	FO0230	FO0231	FO0232	FO0233
FO0234	FO0237	FO0238	FO0239	FO0240	FO0241	FO0242	FO0244	FO0247	FO0249
FO0251	FO0255	FO0257	FO0259	FO0263	FO0264	FO0267	FO0269	FO0271"""
cols = re.split(r'\s+', header_line.strip())
cols.insert(0, "cluster")
pat = re.compile(r'^(.*?)\s*\(\s*([0-9]+(?:\.[0-9]+)?)\s*%\s*\)\s*$')
to_remove = set()
with open(remfile) as rf:
    for line in rf:
        parts = line.strip().split()
        for val in parts[1:]:
            to_remove.add(f"cluster{val}")

print(f"Flagged for removal: {len(to_remove)} clusters")
removed_count = 0

with open(infile, newline='') as fin, open(outfile, 'w', newline='') as fout:
    reader = csv.reader(fin, delimiter='\t')
    writer = csv.writer(fout, delimiter='\t')
    writer.writerow(cols)

    for idx, row in enumerate(reader, start=1):
        label = f"cluster{idx}"
        if label in to_remove:
            removed_count += 1
            continue

        newrow = []
        for cell in row[:len(cols)-1]:
            cell = cell.strip()
            if not cell:
                newrow.append("")
                continue
            m = pat.match(cell)
            if m:
                base, pct = m.group(1), float(m.group(2))
                newrow.append(base if pct >= 95.0 else "")
            else:
                newrow.append(cell)
        if len(newrow) < len(cols)-1:
            newrow.extend([""] * ((len(cols)-1) - len(newrow)))
        newrow.insert(0, label)
        writer.writerow(newrow)
print(f"Removed {removed_count} rows from the output.")
