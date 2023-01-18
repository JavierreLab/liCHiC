import os
import pandas as pd
import sys

## ARGUMENTS
path_vep = sys.argv[1]
path_out = sys.argv[2]
name_general_output = sys.argv[3]
# path_vep = "/home/plopez/Data/GWAS/coding_SNPs_VEP"
# path_out = "/home/plopez/Data/GWAS/VEP_formatted"
# name_general_output = "EUR.phase3.GRCh38.cSNPs_VEP.f.txt"


for filename in os.listdir(path_vep):
    if filename.endswith('VEP.txt'):
        print("Processing file {}...".format(filename))
        vep_df = pd.read_csv(os.path.join(path_vep,filename),sep="\t", skiprows = [x for x in range(0,40)])
        vep_df = vep_df.loc[:,['Location','Gene']]
        vep_df[['chr', 'other']] = vep_df['Location'].str.split(':', 1, expand=True) 
        vep_df[['pos', 'End']] = vep_df['other'].str.split('-', 1, expand=True)
        vep_df.drop('Location', inplace=True, axis=1)
        vep_df.drop('other', inplace=True, axis=1)
        vep_df.drop('End', inplace=True, axis=1) 
        vep_df = vep_df[['chr','pos', 'Gene']]
        vep_df.columns = ['chr','pos', 'ensg']

        ## STORE OUTPUT
        out_name = filename.replace(".txt","_f.txt")
        vep_df.to_csv(os.path.join(path_out,out_name),sep="\t",header=False, index=False)

print("\n")
fout = open(os.path.join(path_out,name_general_output),"w")
fout.write("chr\tpos\tensg\n")

for filename in os.listdir(path_out):
    if filename == name_general_output:
        continue
    print("Joining outputfile {}...".format(filename))
    fd = open(os.path.join(path_out,filename),"r")
    for line in fd:
        fout.write(line)
    fd.close()

fout.close()

