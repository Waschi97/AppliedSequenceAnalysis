from pathlib import Path

scores = str(snakemake.input.scores)
refs = str(snakemake.input.refs)
out = str(snakemake.output)

score_map = {}

with open(scores) as sf:
    for line in sf:
        values = line.split(',')
        genome = str(values[1])
        score = float(values[2])
        if genome in score_map:
            score_map[genome] += [score]
            continue
        score_map[genome] = [score]

best_reference = ""
top_score = float("inf")
for genome, scores in score_map.items():
    avg_score = sum(scores) / len(scores)
    if avg_score >= top_score:
        continue
    top_score = avg_score
    best_reference = genome

#print(f"Best reference: {genome}\navg. blast score: {top_score}")

write = False
for line in open(refs, 'r'):
    if line[0] == ">":
        if write:
            break
        if best_reference.split('|')[1] in line:
            write = True
            Path(out).parent.mkdir(parents=True, exist_ok=True)
            Path(out).touch(exist_ok=False)
    if write:
        with open(out, 'a') as of:
            of.write(line)
