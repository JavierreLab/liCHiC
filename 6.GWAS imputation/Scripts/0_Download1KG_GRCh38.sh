# Download from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/

# download metadata of individuals 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# keep only EUROPEAN individuals
grep EUR integrated_call_samples_v3.20130502.ALL.panel > 1000g_europeans.subjects

# Subset, keep only EUR variants
for i in {1..22}; do
 echo "ALL.chr"$i".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
 vcf-subset -c 1000g_europeans.subjects ALL.chr"$i".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | fill-an-ac | gzip -c > EUR.chr"$i".phase3.GRCh38.vcf.gz
done

# index subsetted VCF files 
gunzip EUR.chr*.vcf.gz
bgzip EUR.chr*.vcf
tabix -p vcf EUR.chr*.vcf.gz
