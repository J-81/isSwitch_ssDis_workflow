import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

#########################
localrules: all, download_scop_ids

configfile: "config/bondugula_V1.yml"



##########################



rule download_pdb_seqs:
    input:
        FTP.remote("ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz")
    output:
        pdb_seqs = "data/external/rcsb/pdb_seqres.txt"
    shell:
        "gunzip -c {input} > {output}"

rule filter_non_proteins:
    """
    Remove all entries except protein entries (contain "mol:protein")
    """
    input:
        rules.download_pdb_seqs.output.pdb_seqs
    output:
        protein_seqs = "data/rcsb.fasta/proteins_seqres.txt"
    log: "logs/non-protein-filter.txt"
    conda: "envs/biopython.yml"
    notebook:
        "notebooks/filter_non_proteins.py.ipynb"

###################

def _scop_url(_version):
    """ get scop parseable file url
    Versions above 1.75 have slighlty different url format
    """
    if _version > 1.75:
        return HTTP.remote(expand("https://scop.berkeley.edu/downloads/parse/dir.des.scope.{_version}-stable.txt", 
                                  _version=_version
                                  ),
                            keep_local=True)
    else:
        return HTTP.remote(expand("https://scop.berkeley.edu/downloads/parse/dir.des.scop.{_version}.txt", 
                                  _version=_version
                                  ),
                            keep_local=True)


rule download_scop_ids:
    input:
        _scop_url(_version=config["scop_version"])
    output:
        expand("data/external/scop.full/dir.des.scop.{_version}.txt", _version=config["scop_version"])
    shell:
        """
        mkdir --parents `dirname {output}` && mv {input} {output}
        """

#####################

rule convert_scop_ids:
    """ Convert to pdb_seqres_style, xxxx_x
    """
    input: 
        rules.download_scop_ids.output
    output:
        expand("data/external/scop.full/scop_{_version}.json", _version=config["scop_version"])
    conda: "envs/biopython.yml"
    notebook: 
        "notebooks/_convert_scop_ids.py.ipynb"

# This is surprisingly fast (completes in minutes on laptop)
rule cluster_pdb:
    input:
        rules.filter_non_proteins.output.protein_seqs
    params:
        prefix=lambda wildcards, output: output.all_seqs[:-15] ,
        min_seq_id = config["min_seq_id"],
    output:
        all_seqs = expand("data/mmseqs2.cluster/cluster{min_seq_id}_all_seqs.fasta", min_seq_id=config["min_seq_id"])[0],
        _ids = expand("data/mmseqs2.cluster/cluster{min_seq_id}_cluster.tsv", min_seq_id=config["min_seq_id"])[0]
    conda: "envs/mmseqs2.yml"
    shell:
        "mkdir -p {params.prefix} && mmseqs easy-cluster {input} {params.prefix} tmp --min-seq-id {params.min_seq_id}"

rule derive_cluster_from_scop:
    """ Find all mmseq2 clusters that include scopPdb5IDs
    """
    input:
        scopPdbIDs = rules.convert_scop_ids.output,
        mmseqsPdbIDs = rules.cluster_pdb.output._ids
    output:
        scop_rep_clusters = expand("data/mmseqs.scop/scop_{scop_version}_clusters.json", scop_version=config["scop_version"])
    conda: "envs/biopython.yml"
    notebook:
        "notebooks/_derive_cluster_from_scop.py.ipynb"




rule extract_clusters:
    """
    From the full cluster file, extract a single multi-sequence fasta
        for each cluster, also cleans up remaining unneeded cluster files
    """
    input:
        all_seqs = rules.cluster_pdb.output.all_seqs # a multi-cluster, multi-sequence fasta
    output:
        fasta = directory("data/mmseqs2.cluster.single")
    conda: "envs/biopython.yml"
    log: "logs/extract_cluster.txt"
    notebook:
        "notebooks/divide_into_individual_fasta_files_per_cluster.py.ipynb"


rule cluster_msa:
    """
    Generate an msa for cluster
        Note: will simply copy fasta file if only single record is present in cluster
    """
    input:
        cluster_fasta = "data/mmseqs2.cluster.single/{_hash}/{clusterID}.fasta"
    output:
        cluster_msa = "data/clustalo.cluster.msa/{_hash}/{clusterID}.fasta",
    conda: "envs/clustalo.yml"
    resources: cpus=28
    shell:
        """
        count=`grep -c ">" {input.cluster_fasta}`; if [ $count -gt 1 ]; then clustalo --threads={resources.cpus} -i {input.cluster_fasta} -o {output.cluster_msa}; else mv {input.cluster_fasta} {output.cluster_msa}; fi
        """

rule download_ssdis:
    input:
        HTTP.remote("https://cdn.rcsb.org/etl/kabschSander/ss_dis.txt.gz")
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output.raw_ssdis)
    output:
        raw_ssdis = "data/external/rcsb/ss_dis.txt"
    shell:
        "mkdir -p {params.output_dir} && gunzip -c {input} > {output}"

rule process_ssdis:
    """ 
    Covert ssdis into a more usable format
    """
    input:
        ss_dis = rules.download_ssdis.output.raw_ssdis
    output:
        reformatted_ss_dis = "data/ss_dis/ss_dis.json"
    conda: "envs/biopython.yml"
    notebook:
        "notebooks/ss_dis_convert.py.ipynb"


rule ssdis_overlay:
    input:
        reformatted_ss_dis="data/ss_dis/ss_dis.json",
        cluster_msa = rules.cluster_msa.output.cluster_msa
    output:
        ssdis_overlay =  "data/cluster.msa.secstruct/{_hash}/{clusterID}.fasta"
    conda: "envs/biopython.yml"
    threads: 4
    notebook:
        "notebooks/overlay_ssdis_onto_msa.py.ipynb"

rule alignment_csv:
    input:
        rules.ssdis_overlay.output.ssdis_overlay
    output:
        ssdis_csv = "data/cluster.csv.secstruct/{_hash}/{clusterID}.csv"
    threads: 4
    conda: "envs/biopython.yml"
    notebook:
        "notebooks/output_alignment_csv.py.ipynb"

rule assess_isSwitch:
    """ Note, this will seek the clusterID file as the output; however,
    it will generate a dataset file for EVERY member of the cluster
    """
    input:
        primary_aln = rules.cluster_msa.output.cluster_msa,
        ssdis_csv = rules.alignment_csv.output.ssdis_csv
    output:
        isSwitchPlus_dataset = "data/isSwitch.dataset/{_hash}/{clusterID}.csv",
        isSwitchPlus_full = "data/isSwitch.full/{_hash}/{clusterID}.csv",
    threads: 4
    conda: "envs/biopython_pandas.yml"
    notebook:
        "notebooks/isSwitch_assessment_plus.py.ipynb"








######################## Control for failed clusters

import json

def _json_to_list(json_path):
    """ A list of mmseqs2 ClusterIDs that contain scop members """
    with open(json_path, "r") as f:
        return json.load(f)

def _hash_clusterIDs(json_path):
    _hash = list()
    clusterIDs = list()
    for item in _json_to_list(json_path):
        _hash.append(item[0:2])
        clusterIDs.append(item)
    return _hash, clusterIDs



def filter_by_dropped(_hash_list, clusterID_list, log_path="logs/dropped.json"):
    with open(log_path, "r") as f:
        failed = json.load(f)

    
    _starting = list(zip(_hash_list, clusterID_list))
    
    filter_out = list()
    for reason, failed_clusterIDs in failed.items():
        # generate hashes associated with ids to filter out
        filter_out_hashes = [_id[0:2] for _id in failed_clusterIDs]

        _remove = list(zip(filter_out_hashes, failed_clusterIDs))
    
        print(f"From {len(_starting)} initial ClusterIDs: ")
        # remove known failed IDs
        for target in _remove:
            _starting.remove(target)
        print(f"\tRemoved {len(_remove)} clusterIDs due to known Reason: {reason}")
    
    surviving_hash, surviving_clusterID = zip(*_starting)
    print(f"After removing known failures: {len(surviving_clusterID)} clusterIDs reached this point")
    return surviving_hash, surviving_clusterID
"""
def _prep_all():
    try:
        _hash, clusterIDs = _hash_clusterIDs(rules.derive_cluster_from_scop.output.scop_rep_clusters[0])
    except:
        print("IF YOU GET THIS ERROR: Run against target rule derive_cluster_from_scop first")
    ############### remove failures ########################
    ## This can be disabled to validate which clusterIDs fail and where
    _hash, clusterIDs = filter_by_dropped(_hash, clusterIDs)
    ############### end remove failures ###################
    # DEBUG
    input("CONFIRM THAT DEBUG IS ON")
    limiter = 10
    _hash = _hash[0:limiter]
    clusterIDs = clusterIDs[0:limiter]
    return expand("data/dataset.final/bondugula_v{version}.csv", version=config["Bondugula_version"])
"""
def _input_aggregate_isSwitch(Wildcards):
    try:
        _hash, clusterIDs = _hash_clusterIDs(rules.derive_cluster_from_scop.output.scop_rep_clusters[0])
    except:
        print("IF YOU GET THIS ERROR: Run against target rule derive_cluster_from_scop first")
    ############### remove failures ########################
    ## This can be disabled to validate which clusterIDs fail and where
    _hash, clusterIDs = filter_by_dropped(_hash, clusterIDs)
    ############### end remove failures ###################
    # DEBUG
    input("CONFIRM THAT DEBUG IS ON")
    limiter = 10
    _hash = _hash[0:limiter]
    clusterIDs = clusterIDs[0:limiter]
    return expand("data/isSwitch.dataset/{_hash}/{clusterID}.csv", zip, _hash=_hash, clusterID=clusterIDs) 


rule aggregate_isSwitch:
    input:
        scop_ids = rules.download_scop_ids.output,
        dataset_csv = _input_aggregate_isSwitch
    output:
        dataset = "data/dataset.final/bondugula_v{version}.csv",
    threads: 4
    conda: "envs/biopython_pandas.yml"
    notebook:
        "notebooks/aggregate_results.py.ipynb"

rule all:
    input:
        expand("data/dataset.final/bondugula_v{version}.csv", version=config["Bondugula_version"])
