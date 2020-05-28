#### Prepare data for Partition Finder
# From a fasta file with multiple samples, this script separates the fasta file into individual fasta files,
# makes a list of sample IDs, separate genes based on a bed file, aligned the genes, merged genes in one file, 
# make a phylip file for partition finder
### You need
# 1. input fasta file
# 2. bedfile "Mtg_template_mod.bed" 
# 3. script: "Fasta2Phylip.pl"

#############

# create bed files for the different genes, tRNA, etc
# based on https://www.ncbi.nlm.nih.gov/nucleotide/NC_005279.1
# bed file starts with 0 - so take 1 out from the start
# bed file is Mtg_template.bed


# To run the script
# 1. Name of script
# 2. Name of input fasta file
# 3. Number of sequences in your fasta file
# 4. Number of genes in the mitogenome. Usually 38
# 5. Directory path
# 6. Name of output file
# Example: ./Fasta2aligned_genes.sh Bmys-ancmod_30samples16393bp_svalbard 30 38 /shared/volume/hologenomics/data/andreaac/bowhead/partition_finder Bmys_ancmod_30ind_aligned

datainfile=$1
n_sequences=$2
n_genes=$3
directory_path=$4
alingned_filename=$5


#separate one fasta into multiple fasta - one fasta for each individual
cat ${directory_path}/${datainfile}.fasta | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")}
        print $0 > filename
}'

# make a list of all sample ids
awk 'sub(/^>/, "")' ${datainfile}.fasta > samplelist.txt


# create a new folder
mkdir separate_genes
# go in that folder
cd ${directory_path}/separate_genes


#make y number of copies of each line, the y represents the different genes from the mitogenome (usually 38)
awk '{for(i=0;i<38;i++)print}' ${directory_path}/samplelist.txt > samplelist_${n_genes}times.txt

# make x (number of sequences) number of copies of bed file after each other without column 1
for i in $(seq 1 ${n_sequences}); do cat ${directory_path}/Mtg_template_mod.bed >> Mtg_template_${n_sequences}times.bed; done

#remove column 1
awk ' {print$2"\t"$3"\t"$4"\t"$5"\t"$6} ' Mtg_template_${n_sequences}times.bed > Mtg_template_${n_sequences}times_column2-6.bed

#paste the sample id column with long bed file.
paste samplelist_${n_genes}times.txt Mtg_template_${n_sequences}times_column2-6.bed > Mtg_${n_sequences}_01.bed


mv ${directory_path}/*.fasta ${directory_path}/separate_genes
mv ${directory_path}/separate_genes/${datainfile}.fasta ${directory_path}/


##### 

#load dos2unix
module load dos2unix/v7.3.4

#remove windows caracters
dos2unix *.fasta

# use bedtools to seperate genes
#load bedtools
module load bedtools/v2.26.0

#Separate genes using bedtools

cd ${directory_path}/separate_genes
for file in *.fasta
do bedtools getfasta -name -s -fi ${file} -bed Mtg_${n_sequences}_01.bed > ${file}_SepGenes.fa
done

# combine genes from all individuals
#You need to add the name of only one of the files in order to get a list of genes

mkdir ${directory_path}/separate_genes/unaligned_genes

for gene in `grep ">" Bmys_17-19.fasta_SepGenes.fa | cut -d ":" -f 1 | sed ' s/>//g ' `
do
for file in *.fasta_SepGenes.fa
do grep ">${gene}:" -A1 $file >> unaligned_genes/${gene}.fasta
done
done


#load mafft
module load mafft/v7.310

cd ${directory_path}/separate_genes/
mkdir ${directory_path}/separate_genes/aligned_genes



#Align genes
cd ${directory_path}/separate_genes/unaligned_genes

for file in *fasta
do
xsbatch -c 1 --time 01:00:00 --mem-per-cpu 5000 -- "mafft --maxiterate 1000 --globalpair $file > ../aligned_genes/${file%.*}.aln.fasta"
done

sleep 60

#Join genes for all regions

#remove windows caracters
dos2unix ${directory_path}/separate_genes/aligned_genes/*.fasta

cd ${directory_path}/separate_genes/aligned_genes
rm ${alingned_filename}.fasta

for sample in `grep ">" CDS_COX1.aln.fasta | cut -f3 -d":"`
do
echo ">${sample}" >> ${alingned_filename}.fasta
for gene in ATP6 ATP8 COX1 COX2 COX3 CYTB ND1 ND2 ND3 ND4 ND4L ND5 ND6 Ctrl-Reg_D-loop 12S 16S Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu_1 Leu_2 Lys Met Phe Pro Ser_1 Ser_2 Thr Trp Tyr Val
do
perl -pe 'if ($_=~/^>/){s/\n/\t/};s/\n//g;s/>/\n>/g;s/\t/\n/g' *$gene.aln.fasta | tail -n+2 | grep "$sample" -A1 | tail -n+2 >> ${alingned_filename}.fasta
done
done


#make phylip alignment
#perl  Fasta2Phylip.pl ${alingned_filename}.fasta ${alingned_filename}.phylip

perl ${directory_path}/Fasta2Phylip.pl ${directory_path}/separate_genes/aligned_genes/${alingned_filename}.fasta ${directory_path}/separate_genes/aligned_genes/${alingned_filename}.phylip

#make directory for partition analysis

mkdir ${directory_path}/partition_finder_analysis

cp ${directory_path}/separate_genes/aligned_genes/${alingned_filename}.phylip ${directory_path}/partition_finder_analysis/
