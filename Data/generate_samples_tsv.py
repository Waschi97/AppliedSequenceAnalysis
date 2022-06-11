import os
from pathlib import Path


fastq_directory = "/storage/mi/tomw97/Data/AppliedSeqAn/contaminants/reads_with_contaminants/fastq"
db_directory = "/storage/mi/tomw97/Data/AppliedSeqAn/contaminants/reads_with_contaminants/info"
out_file = "/storage/mi/tomw97/Development/AppliedSequenceAnalysis/Data/project2_samples.tsv"
# tailing chars after the sample name (without the file ending(s))
x = 12

with open(out_file, 'w') as of:
    of.write("sample\tfq1\tfq2\tc_db\n")
    for file in os.listdir(fastq_directory):
        # skip non fastq.gz files
        if not file.endswith(".fastq.gz"):
            continue

        # remove file endings
        file_name = Path(file)
        while file_name.suffix != '':
            file_name = file_name.with_suffix('')

        if str(file_name)[-1] == '1':
            of.write(f"{str(file_name.stem)[:-x]}\t{fastq_directory}/{file}\t{fastq_directory}/{str(file_name)[:-1]}2.fastq.gz\t")

            db_found = False
            for db in os.listdir(db_directory):
                if db.endswith('.fasta') and str(file_name.stem)[:-x] in str(db):
                        of.write(f"{db_directory}/{db}\n")
                        db_found = True
                        break
            if not db_found:
                of.write(f"-\n")
            

