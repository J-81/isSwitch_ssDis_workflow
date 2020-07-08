""" Reads mmseqs2 style all_seqs fastsa and generates an MSA on all clusters with > 1 sequence"""

import os
import subprocess

MIN_CLUSTER_SIZE = 2
MAX_CLUSTER_SIZE = 5000

def parse_validate_new_cluster_id(line):
    new_cluster = line.replace(">","").strip()
    if len(new_cluster) < 6:
        return None
    else:
        return new_cluster[0:4] + new_cluster[4] + new_cluster[5].lower() # ugly way to lowercase last letter

# from https://www.geeksforgeeks.org/python-check-two-lists-least-one-element-common/
def common_member(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if len(a_set.intersection(b_set)) > 0: 
        return list(a_set.intersection(b_set))[0]
    return None

os.makedirs(snakemake.output.cluster_msa, exist_ok=True)
cluster_buffer = ""
new_cluster = False
cur_cluster_ids = list()

# parse bondugla ids to xxxx_x format
# original format example: d1ejga_
with open(snakemake.input.bondugula, "r") as f:
    bondugula_ids = [f"{_id[1:5]}_{_id[5]}"  for _id in f.readlines()]

with open(snakemake.input.all_fasta, "r") as f:
    for line in f.readlines():
        ##### NEW CLUSTER LINE ######
        if len(line) < 10: # all full desc are at least this length

            if not new_cluster:
                new_cluster = parse_validate_new_cluster_id(line)
                continue

            outpath = os.path.join(snakemake.output.cluster_msa,new_cluster[0:2],f"{new_cluster}.fasta")
            

            ######### WRITE AND MSA ##############
            if all((MAX_CLUSTER_SIZE >= cluster_buffer.count(">") >= MIN_CLUSTER_SIZE,
                    common_member(bondugula_ids, cur_cluster_ids))): # skip any clusters of size 1
                bond_id = common_member(bondugula_ids, cur_cluster_ids)
                os.makedirs(os.path.dirname(outpath), exist_ok=True)
                print(new_cluster, bond_id)
                bondugula_ids.remove(bond_id)
                tmp = "tmp/clustalo.fasta"
                with open(tmp, "w") as f:
                    f.write(cluster_buffer)
                subprocess.run(f"clustalo --threads={snakemake.threads} -i {tmp} -o {outpath}", shell=True)
            
            
            ######## RESET BUFFER
            cluster_buffer = ""
            new_cluster = parse_validate_new_cluster_id(line)
            cur_cluster_ids = list()

            continue

        ##### NON-NEW  CLUSTER LINE
        cluster_buffer += line.replace("'","") # remove all single quotes, these will break the subprocess command
        if ">" in line:
            cur_cluster_ids.append(parse_validate_new_cluster_id(line))

os.makedirs(os.path.dirname(snakemake.output.missing), exist_ok=True)
with open(snakemake.output.missing, "w") as f:
    [f.write(_id+"\n") for _id in bondugula_ids]
