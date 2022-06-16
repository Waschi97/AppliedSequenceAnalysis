from pathlib import Path
from shutil import copyfile

score_map = {}

with open(snakemake.input.scores) as sf:
    for line in sf:
        values = line.split('\t')
        genome = str(values[1])
        score = float(values[2])
        if genome in score_map:
            score_map[genome] += [score]
            continue
        score_map[genome] = [score]

best_reference = ""
top_score = float("-inf")
for genome, scores in score_map.items():
    avg_score = sum(scores) / len(scores)
    if avg_score <= top_score:
        continue
    top_score = avg_score
    best_reference = genome

#print(f"Best reference: {genome}\navg. blast score: {top_score}")

for ref_file in snakemake.input.refs:
    if best_reference in ref_file:
        best_reference = ref_file

copyfile(best_reference, str(snakemake.output))
