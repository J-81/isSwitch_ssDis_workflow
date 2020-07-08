import os

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()




rule download_pdb_seqs:
    input:
        FTP.remote("ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz")
    output:
        "data/external/pdb_seqres.txt"
    shell:
        "gunzip -c {input} > {output}"

# NOTE: 25% identity cutoff used hardcoded in url here
# this can be parameterized if other options are wanted
rule download_bondugula_ids:
    output:
        "data/external/bondugula_ids.txt"
    shell:
        "wget -O {output} 'https://scop.berkeley.edu/astral/subsets/?ver=1.75&get=bib&seqOption=0&item=ids&cut=25'"

rule download_ssdis:
    input:
        HTTP.remote("https://cdn.rcsb.org/etl/kabschSander/ss_dis.txt.gz")
    output:
        "data/external/ss_dis.txt"
    shell:
        "gunzip -c {input} > {output}"

rule process_ssdis:
    input:
        ss_dis="data/external/ss_dis.txt"
    output:
        reformatted_ss_dis="data/ss_dis/ss_dis.json"
    conda: "envs/biopython.yml"
    notebook:
        "notebooks/ss_dis_convert.py.ipynb"

rule cluster_pdb:
    input:
        "data/external/pdb_seqres.txt"
    params:
        prefix="data/mmseqs2.cluster/cluster98"
    output:
        "data/mmseqs2.cluster/cluster98_all_seqs.fasta"
    conda: "envs/mmseqs2.yml"
    shell:
        "mmseqs easy-cluster {input} {params.prefix} tmp --min-seq-id 0.98"

rule cluster_msa:
    input:
        all_fasta="data/mmseqs2.cluster/cluster98_all_seqs.fasta",
        bondugula="data/external/bondugula_ids.txt"
    output:
        cluster_msa = directory("data/cluster.msa"),
        missing = "data/filtered/missing_from_all_seqs.txt"
    conda: "envs/clustalo.yml"
    threads: 4
    script:
        "scripts/msa_from_cluster_fasta.py"



_hash, _pdb5id = glob_wildcards("data/cluster.msa/{hash}/{pdb5ID}.fasta")
rule ssdis_overlay:
    input:
        reformatted_ss_dis="data/ss_dis/ss_dis.json",
        cluster_msa = expand("data/cluster.msa/{hash}/{pdb5ID}.fasta", zip, hash=_hash, pdb5ID=_pdb5id)
    output:
        ssdis_overlay =  expand("data/cluster.msa.secstruct/{hash}/{pdb5ID}.fasta", zip, hash=_hash, pdb5ID=_pdb5id),
        missing_ssdis = "data/filtered/MISSING_SSDIS.txt"
    conda: "envs/biopython.yml"
    threads: 4
    notebook:
        "notebooks/overlay_ssdis_onto_msa.py.ipynb"

_hash, _pdb5id = glob_wildcards("data/cluster.msa.secstruct/{hash}/{pdb5ID}.fasta")
rule alignment_csv:
    input:
        rules.ssdis_overlay.output.ssdis_overlay
    output:
        ssdis_csv = expand("data/cluster.csv/{hash}/{pdb5ID}.csv", zip, hash=_hash, pdb5ID=_pdb5id)
    threads: 4
    conda: "envs/biopython.yml"
    notebook:
        "notebooks/output_alignment_csv.py.ipynb"

_hash, _pdb5id = glob_wildcards("data/cluster.csv/{hash}/{pdb5ID}.csv")
rule assess_isSwitch:
    input:
        primary_aln = rules.ssdis_overlay.input.cluster_msa,
        ssdis_csv = rules.alignment_csv.output.ssdis_csv
    output:
        isSwitchPlus_dataset = directory("data/isSwitch.dataset"),
        isSwitchPlus_full = expand("data/isSwitch.full/{hash}/{pdb5ID}.csv", zip, hash=_hash, pdb5ID=_pdb5id),
        missing = "data/filtered/missing_due_to_missing_ssdis.txt"
    threads: 4
    conda: "envs/biopython_pandas.yml"
    notebook:
        "notebooks/isSwitch_assessment_plus.py.ipynb"

rule aggregate_isSwitch:
    input:
        bondugula_ids = "data/external/bondugula_ids.txt",
        cluster_csv = rules.assess_isSwitch.output.isSwitchPlus_dataset
    output:
        dataset = "data/isSwitch.full.aggregate/bondugula.csv",
        missing = "data/filtered/missing_in_final_aggregation.txt"
    threads: 4
    conda: "envs/biopython_pandas.yml"
    notebook:
        "notebooks/aggregate_results.py.ipynb"


# a target rule to define final output goal
rule all:
    input:
        rules.aggregate_isSwitch.output.dataset
