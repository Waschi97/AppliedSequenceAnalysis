from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import math

def SeqVar(basecount, N):
    variance = 0
    for c in basecount:
        if c == 0:
            continue
        variance += float(c)/N * math.log(float(c)/N)
    return variance * -1

window_size = snakemake.params.window_size

is_initialised = False

for record in SeqIO.parse(snakemake.input[0], "fasta"):

    if not is_initialised:
        # initialise the dataframe with zeros
        df = np.zeros((len(record.seq), 4),dtype=np.uint64)
        is_initialised = True
    else:
        # check if a sequence has another length (i.e. is not aligned)
        if len(record.seq) != df.shape[0]:
            raise Exception("Not all sequences are of the same length! Recalculate your alignment!")

    # count occurences
    for i in range(0,len(record.seq)):
        base = str(record.seq[i]).lower()
        if base == 'a':
            df[i, 0] += 1
            continue
        if base == 'c':
            df[i, 1] += 1
            continue
        if base == 'g':
            df[i, 2] += 1
            continue
        if base == 't':
            df[i, 3] += 1
            continue

# calculate sequence variance for each position
seq_var = [SeqVar(row,df.shape[0]) for row in df]

# apply window size
num_windows = math.floor(len(seq_var)/window_size)
seq_var_windowed = [sum(seq_var[window_size*i:window_size*(i+1)])/window_size for i in range(0, num_windows)]

# handle last window which might be smaller than window_size
remainder = len(seq_var) % window_size
if remainder != 0:
    seq_var_windowed = seq_var_windowed + [sum(seq_var[-remainder:])/remainder]

# write txt
X = [] # save window start for plotting
with open(snakemake.output.txt, 'w') as f:
    for i in range(0, len(seq_var_windowed) -1):
        start = i * window_size + 1
        end = (i+1) * window_size
        f.write(f"{start}-{end}: {seq_var_windowed[i]}\n")
        X.append(start)

    start_last_window = window_size * num_windows + 1
    f.write(f"{start_last_window}-{len(seq_var)}: {seq_var_windowed[-1]}")
    X.append(start_last_window)

# plotting
fig = plt.figure()
ax = plt.axes()
plt.plot(X, seq_var_windowed)
ax.set_xlabel("window starting position")
ax.set_ylabel("avg. sequence variability")
ax.set_title(f"Sequence Variability for Windows of Size {window_size} bp")
fig.tight_layout()
fig.savefig(snakemake.output.png)