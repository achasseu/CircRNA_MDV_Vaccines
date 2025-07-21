####python script to obtain the stranding of the reads in the sam file
import pandas as pd

split_sam_file = "../../Downloads/OneDrive_2025-07-16/Analyse_Camille/GaHV3/aln_split.sam" #### adapt this line for your own file
fasta_file = "../../Downloads/OneDrive_2025-07-16/Analyse_Camille/GaHV3/GaHV3_SB1.fasta" #### adapt this line for your own file
#read the sam file
with open(split_sam_file, 'r') as sam_file:
    sam_file=sam_file.readlines()
with open(fasta_file, 'r') as fasta_file:
    fasta_file=fasta_file.readlines()
    fasta_file=fasta_file[1:]
    fasta_file="".join(fasta_file)

### function to extract the real sequence from the mapping
def extract_sequence(mappings):
    df=pd.DataFrame(mappings)
    df.columns=["read_id","flag","chromosome","position","mapq","cigar","mate_chromosome","mate_position","other unknown","sequence","quality","insert_size","template_length","number_of_mappings","mapping_quality","mapping_name"]    
    if "0" in df["flag"].values or "16" in df["flag"].values:
        df=df.sort_values(by="position")
        df=df.sort_values(by="flag")
        sequence_totale=df[df["flag"].isin(["0", "16"])]["sequence"].values[0]
        df=df.values.tolist()
        for mapping in df[1:]:
            sequence_totale=sequence_totale.replace(mapping[9],"")
        df[0][9]=sequence_totale
        df=pd.DataFrame(df)
        df.columns=["read_id","flag","chromosome","position","mapq","cigar","mate_chromosome","mate_position","other unknown","sequence","quality","insert_size","template_length","number_of_mappings","mapping_quality","mapping_name"]
        df=df.sort_values(by="position")
        df=df.values.tolist()
        return df
    else:
        return None
  

### function to extract the splice junctions from the mappings
def extract_splice_junctions(mappings,fasta_file):
    index=0
    sense=0
    df=pd.DataFrame(mappings)
    df.columns=["read_id","flag","chromosome","position","mapq","cigar","mate_chromosome","mate_position","other unknown","sequence","quality","insert_size","template_length","number_of_mappings","mapping_quality","mapping_name"]    
    read_id=df["read_id"].values[0]
    if "0" in df["flag"].values:
        sense=1
    else:
        sense=-1
    for mapping in mappings:
        if mapping!=mappings[-1]:
            coordinates=[int(mapping[3])+len(mapping[9]),int(mappings[index+1][3])]
            if coordinates[0]>99638:
                coordinates[0]=coordinates[0]+101
                coordinates[1]=coordinates[1]+101
            splice_junction=fasta_file[coordinates[0]-5:coordinates[0]+5]+fasta_file[coordinates[1]-5:coordinates[1]+5]
            index+=1
            return [read_id,splice_junction,sense]
        else:
            return None

### determining the sense and the number of splice junctions for each read

##### create lists of reads that have the same read ID (since they are split)
mappings_list=[]
read_ids=[]
index=-1
print("--------------------------------")
print("Starting the script")
print("--------------------------------")
print("First, we will create a list of splice junctions")
print(len(sam_file),"lines in the sam file",sep=" ")
splice_junctions_list=[]
for line in sam_file: ###create a list of lists of reads that have the same read ID
    mappings=[]
    index+=1
    if line.startswith('@'):
        continue
    else:
        read_id=line.split("\t")[0]
    
    if not read_id in read_ids[-5:]:
        mappings.append(line.replace(";\n","").split("\t"))
        for second_line in sam_file[index+1:index+10]:
            if second_line.replace(";\n","").split("\t")[0]==read_id:
                mappings.append(second_line.replace(";\n","").split("\t"))
    else:
        continue
    if len(mappings)>1:
        mappings=extract_sequence(mappings)
        if mappings is not None:
            splice_junctions_list.append(extract_splice_junctions(mappings,fasta_file))

    read_ids.append(read_id)
    if index%10000==0:
        print(round((index/len(sam_file))*100, 2))

splice_junctions_df = pd.DataFrame(splice_junctions_list, columns=["read_id", "splice_junction", "sense"])

print("--------------------------------")
print("Done")
print("--------------------------------")
print("Now, we will count the splice junctions and orient the misstranded reads")

# Group by splice_junction, aggregate read_ids into a list, and calculate average sense. 
# Store the length of the read_ids list in a new column called "count". Sort by count.
result_df = (
    splice_junctions_df
    .groupby("splice_junction")
    .agg(
        read_ids=("read_id", lambda x: list(x)),
        avg_sense=("sense", "mean"),
        count=("read_id", "size")
    )
    .reset_index()
    .sort_values(by="count", ascending=False)
)

#if avg_sense is between -0.5 and 0.5, we'll check the splice_junction to try to determine the sense
def determine_sense(splice_junction,avg_sense):

    if "GT" in splice_junction[:(len(splice_junction)//2)+1] and "AG" in splice_junction[len(splice_junction)//2:]:
        return 2
    elif "CT" in splice_junction[:(len(splice_junction)//2)+1] and "AC" in splice_junction[len(splice_junction)//2:]:
        return -2
    else:
        return avg_sense

result_df["avg_sense"]=result_df.apply(lambda row: determine_sense(row["splice_junction"],row["avg_sense"]), axis=1)

#Now, we'll extract all the read_ids from the rows that have a count value superior to 2
#Then, we'll extract the sequences from the sam file into alignment_sense.sam and alignment_antisense.sam for positive and negative sense respectively
print("--------------------------------")
print("Done")
print("--------------------------------")
print("Now, we will filter the reads based on their orientation and their count >2")

read_ids_sense=result_df[result_df["count"]>2]
read_ids_sense=read_ids_sense[read_ids_sense["avg_sense"]>0]["read_ids"].values
read_ids_sense = [read_id for sublist in read_ids_sense for read_id in sublist]
read_ids_antisense=result_df[result_df["count"]>2]
read_ids_antisense=read_ids_antisense[read_ids_antisense["avg_sense"]<0]["read_ids"].values
read_ids_antisense = [read_id for sublist in read_ids_antisense for read_id in sublist]


print("--------------------------------")
print("Done")
print("--------------------------------")
print("Now, we will create the alignment files")

sense_ids = set(read_ids_sense)         
antisense_ids = set(read_ids_antisense) 

sense_file=open("../../Downloads/OneDrive_2025-07-16/Analyse_Camille/GaHV3/aln_sense.sam", "w") #### adapt this line for your own file
antisense_file=open("../../Downloads/OneDrive_2025-07-16/Analyse_Camille/GaHV3/aln_antisense.sam", "w") #### adapt this line for your own file
index=0
for index, line in enumerate(sam_file):
    readid=line.split("\t")[0]
    if readid in sense_ids:
        sense_file.write(line)
    elif readid in antisense_ids:
        antisense_file.write(line)
    if index%10000==0:
        print(round((index/len(sam_file))*100, 2))
    index+=1
sense_file.close()
antisense_file.close()

print("--------------------------------")
print("Done")
print("--------------------------------")
print("The alignment files have been created")
print("--------------------------------")
print("{} reads have been filtered".format(len(read_ids_sense)+len(read_ids_antisense)))
print("--------------------------------")
print("Well done!")
print("--------------------------------")