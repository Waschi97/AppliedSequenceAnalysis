from Bio import Phylo
import matplotlib.pyplot as plt

tree = Phylo.read(str(snakemake.input), "newick")
f = plt.figure()
Phylo.draw(tree)
plt.savefig(str(snakemake.output))