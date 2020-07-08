from Bio import SeqIO
import json


def secstr_disorder_merge(secstr, disorder):
    #replace all spaces with dashes
    replacement = list(secstr.replace(" ", "-"))
    #if there should be an X, replace the dash with an X
    for i in range(len(replacement)):
        if disorder[i] == "X":
            replacement[i] = "X"
    replacement_str = ''.join(replacement)
    return replacement_str

def create_dict():
    data = {}
    record_list = []

    print("Starting to read file and write into ss_dis-2.txt...")
    #replace all the spaces in the file and replaces it with dashes
    #so that no information goes missing when parsing through the file
    replace = open("tmp/ss_dis_interim.txt", "w+")
    with open(snakemake.input.ss_dis,"r+") as in_file:
        for line in in_file:
            fixed_line = line.replace(" ", "-")
            replace.write(fixed_line)
    replace.close()
    print("Success")

    print("Now to parse through the whole thing....")
    # After all the spaces have been replaced, move on
    _records = [r for r in SeqIO.parse("tmp/ss_dis_interim.txt", "fasta")]

    record_dict = SeqIO.to_dict(_records)
    
    for i, record in enumerate(_records):
        print(f"Fetching record {i+1} of {len(_records)} {' '*20}", end="\r")
        #creates unique id for dictionary
        id = record.id.split(":")[0] + ":" + record.id.split(":")[1]
        #if unique, add to list
        if(not (id in record_list)):
            record_list.append(id)
    print("Success")

    print("Finally, to add everything into a dictionary...")
    for i, id in enumerate(record_list):
        print(f"Converting {i+1} of {len(record_list)} {' '*8}", end="\r")
        #create new id for future use
        new_id = id.split(":")[0] + id.split(":")[1]
        #find the sequence, secstr, and disorder sequences from dictionary
        sequence = str(record_dict[id + ":sequence"].seq)
        secstr = str(record_dict[id + ":secstr"].seq)
        disorder = str(record_dict[id + ":disorder"].seq)
        merge = str(secstr_disorder_merge(secstr, disorder))
        #add all information into useful dictionary
        data[new_id] = sequence, merge
    print("Success")

    print("Finally to create the JSON File...")
    with open(snakemake.output.reformatted_ss_dis, "w+") as file:
            json.dump(data,file, indent=4)
    print("Finished!")



create_dict()
