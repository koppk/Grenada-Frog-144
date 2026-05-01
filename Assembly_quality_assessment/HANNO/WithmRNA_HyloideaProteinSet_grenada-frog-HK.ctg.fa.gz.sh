#!/bin/bash
# Author: Kopp K, Pristimantis euphronides genome project
set -e
#set -o pipefail

export SCRIPTS=/data/software/HANNO/scripts

##RUN HANNO assembly pipeline:

mkdir WithmRNA_HyloideaProteinSet_grenada-frog-HK.ctg.fa.gz
cd WithmRNA_HyloideaProteinSet_grenada-frog-HK.ctg.fa.gz
##MAP PROTEINS in ../../GrenadaFrog144/Hyloidea_Proteins/combined_Hyloidea_proteins.faa
date
seqtk seq ../grenada-frog-HK.ctg.fa.gz > asm.fa
seqtk seq ../../../GrenadaFrog144/Hyloidea_Proteins/combined_Hyloidea_proteins.faa > protein.input.fa
miniprot --outn=10 -p 0.5 --outc=0.5 --outs=0.5 -Iu -t 24 --gtf asm.fa protein.input.fa | gtfToGenePred stdin stdout | genePredToBed stdin stdout | awk -f /data/software/HANNO/scripts/bed12ToGTF_addscore-tacoFake.awk > asm.prot.gtf
##MERGE USING STRINGTIE
stringtie --merge -o prot.merge.gtf asm.prot.gtf
##MERGE USING TACO

ls asm.prot.gtf > TACO.gtf.list

mamba run -n TACO taco_run --max-isoforms 7 -p 24 -o output1 TACO.gtf.list


##GET BESTMODELS FROM CLUSTER OF PREDICTIONS
awk 'BEGIN{OFS="\t";FS="\t"} {if($3=="exon"){$3="CDS";print}}' prot.merge.gtf > prot.merge1.cds.gtf
bash /data/software/HANNO/scripts/BESTCLUSTERGTF-SEQOUT-e.sh prot.merge1.cds.gtf 0.1 asm.fa
awk 'BEGIN{OFS="\t";FS="\t"} {if($3=="exon"){$3="CDS";print}}' output1/assembly.gtf > prot.merge2.cds.gtf
bash /data/software/HANNO/scripts/BESTCLUSTERGTF-SEQOUT-e.sh prot.merge2.cds.gtf 0.1 asm.fa
cat prot.merge1.cds.gtf.clustered.gtf prot.merge2.cds.gtf.clustered.gtf > MODELS1.gtf 
cat prot.merge1.cds.gtf.clustered.bed12 prot.merge2.cds.gtf.clustered.bed12 > MODELS1.bed12
cat prot.merge1.cds.gtf.clustered.cds.faa prot.merge2.cds.gtf.clustered.cds.faa > MODELS1.faa
rm -f prot.merge1* prot.merge2*

##MAP mRNAs in ../../GrenadaFrog144/Hyloidea_Proteins/combined_rna.fna.gz		
seqtk seq ../../../GrenadaFrog144/Hyloidea_Proteins/combined_rna.fna.gz > rna.input.fa
awk 'BEGIN{OFS="\t";FS=OFS;}{$3="exon";print;}' MODELS1.gtf | paftools.js gff2bed -j - > splicejunctions.bed
minimap2 -I 100G -t 24 -x splice -a --junc-bed splicejunctions.bed asm.fa rna.input.fa > rna.input.sam
samtools sort -o rna.input.bam rna.input.sam
rm -f rna.input.sam
##convert to bed12 and correct strands and check for bad transcripts ("exon length 0")
samtools view rna.input.bam | awk 'BEGIN{OFS="\t";FS="\t";} {counter[$1]++;$1=$1"_"counter[$1];n=split($0,a,"\t");printf $1"\t";for(x=1;x<=n;x++){if(substr(a[x],1,3)=="ts:"){printf a[x]}};printf "\n"}' | grep "ts:A:" | sed "s/ts:A://g" > rna.strands
##BUG FIX HERE: bamToBed output misses commas at end of field 11 and 12!!! Needs Correction by awk script, otherwise bed12ToGTF.awk or bed12ToGTF_addscore-tacoFake.awk will miss terminal exons, which leads to missing UTRs!
bamToBed -split -bed12 -i rna.input.bam | awk 'BEGIN{OFS="\t";FS=OFS;}{$11=$11",";$12=$12",";print}' > rna.bed12
awk -v infile=rna.strands 'BEGIN{OFS="\t";FS="\t";while(getline l < infile){split(l,a,"\t");h[a[1]]=a[2]}} {if($10==1){frame="."} else{frame=$6};counter[$4]++;$4=$4"_"counter[$4];if(h[$4]=="-" && h[$4]!=""){if($6=="+"){$6="-"} else{$6="+"}};if(frame=="."){$6="."};print}' rna.bed12 | awk '{n=split($11,a,",");for(x=1;x<n;x++){if(a[x]==0){t=1} else{t=0}};if(t==0){print}}' > rna.stranded.bed12

##merge protein mappings + denovo transcritome using TACO and stringtie (add gtf here if provided)
##TACO

cat MODELS1.bed12 rna.stranded.bed12 | awk -f /data/software/HANNO/scripts/bed12ToGTF_addscore-tacoFake.awk > evidence.stringtie-like.gtf
ls evidence.stringtie-like.gtf  > TACO.gtf.list2

mamba run -n TACO taco_run --max-isoforms 7 -p 24 -o output2 TACO.gtf.list2

##Stringtie
stringtie --merge -o evidence.merge.gtf TACO.gtf.list2
stringtie --merge -G output2/assembly.gtf -o MODELS2.gtf evidence.merge.gtf output2/assembly.gtf
sed "s/MSTRG./MSTRG/g" MODELS2.gtf | awk 'BEGIN{OFS="\t";FS=OFS;}{gsub(/\./,"_",$9);gsub("_p",".p",$9);print}' > MODELS2.gtf.t 
mv MODELS2.gtf.t MODELS2.gtf 
gtfToGenePred MODELS2.gtf stdout | genePredToBed stdin MODELS2.bed12
##Transdecoder ORFs
gtf_to_alignment_gff3.pl MODELS2.gtf > MODELS2.gff
gtf_genome_to_cdna_fasta.pl MODELS2.gtf asm.fa > MODELS2.mrna.fa
seqtk seq -l 0 MODELS2.mrna.fa | split -l 4000 -
ls x?? | awk '{print "TransDecoder.LongOrfs -S -m 80 -t "$1}' | parallel -j 24 > parallel.ORFs.log 2>&1
cat x??.transdecoder_dir/longest_orfs.pep  > MODELS2.longest_orfs.pep
cat x??.transdecoder_dir/longest_orfs.gff3  > MODELS2.longest_orfs.gff3
rm -rf x?? x??.transdecoder_dir

seqtk seq ../tetrapoda_odb10/refseq_db.faa > proteins.reference.faa

##score ORFs by last alignment to references
lastdb -P 24 -p PROTDB proteins.reference.faa
lastal -P 24 -p BL80 -m100 -K 1 PROTDB MODELS2.longest_orfs.pep > MODELS2.maf
maf-convert blasttab MODELS2.maf| sort -k1,1V -k12,12rn| awk '{n=split($1,a,".");name=a[1];for(x=2;x<n;x++){name=name"."a[x]};print name"\t"$0}'| awk '{if(o!=$1){print};o=$1}' > MODELS2.bestORF
touch MODELS2.bestORF.remove MODELS2.bestORF.keeplist
awk '{if($13>=200){print $2 > "MODELS2.bestORF.keeplist";print $1".p" > "MODELS2.bestORF.remove"}}' MODELS2.bestORF
##Create list of best ORFs
awk '{if($3=="CDS") {split($0,d,"=");print $1"\t"$5-$4+1"\t"d[3]}}' MODELS2.longest_orfs.gff3 | sort -k1,1 -k2,2rn | awk '{if(oldid!=$1) {split($3,a,";");print a[1]};oldid=$1}' > MODELS2.transdecoder.keeplist
grep -vFf MODELS2.bestORF.remove MODELS2.transdecoder.keeplist | cat - MODELS2.bestORF.keeplist > MODELS2.transdecoder.keeplist2
grep -w -F -f MODELS2.transdecoder.keeplist2 MODELS2.longest_orfs.gff3 > MODELS2.longest_orfs.gff3_BEST2
##shift CDS to ATG codons
grep -w CDS MODELS2.longest_orfs.gff3_BEST2 | awk '{print $1"\t"$4-1"\t"$5}' > MODELS2.start.bed
bedtools getfasta -split -s -bed MODELS2.start.bed -fi MODELS2.mrna.fa -fo MODELS2.start.fa
awk '{if(substr($1,1,1)==">"){split($0,a,/[>:]/)} else{seq=toupper($1);for(x=1;x<=length(seq);x=x+3){if(substr(seq,x,3)=="ATG"){n=x;break}};print a[2]"\t"n-1; }}' MODELS2.start.fa > MODELS2.start.shift
awk -v infile=MODELS2.start.shift -v maxshiftpercent=20 'BEGIN{OFS="\t";FS="\t";while(getline l < infile){split(l,a,"\t");h[a[1]]=a[2]}} {if($3=="CDS" && 100*h[$1]/($5-$4-1)<=maxshiftpercent &&  h[$1]>0){$4=$4+h[$1];shift++};print} END{print "shifted "shift" start-codons!" > "/dev/stderr"}' MODELS2.longest_orfs.gff3_BEST2 > MODELS2.longest_orfs.gff3_BEST2.corrATG
##transfer ORF coordinates to genome
cdna_alignment_orf_to_genome_orf.pl MODELS2.longest_orfs.gff3_BEST2.corrATG MODELS2.gff MODELS2.mrna.fa > MODELS2.CDS2.gff3
awk -f  /data/software/HANNO/scripts/TRANSDECODER-newversion-GFF3toGTF.awk MODELS2.CDS2.gff3 > MODELS2.CDS2.gtf
gtfToGenePred MODELS2.CDS2.gtf MODELS2.CDS2.gp
genePredToBed MODELS2.CDS2.gp MODELS2.CDS2.bed12
bash /data/software/HANNO/scripts/CDS_gtfToBed12 MODELS2.CDS2.gtf > MODELS2.CDS2only.bed12
bedtools getfasta -split -name -s -bed MODELS2.CDS2.bed12 -fi asm.fa -fo /dev/stdout > MODELS2.CDS2.mRNA.fa
bedtools getfasta -split -name -s -bed MODELS2.CDS2only.bed12 -fi asm.fa -fo /dev/stdout > MODELS2.CDS2.fa
/data/software/HANNO/scripts/TRANSLATE.sh MODELS2.CDS2.fa > MODELS2.CDS2.faa
##CLEAN-UP
ls MODELS2* | grep -v MODELS2.CDS2 | xargs rm

##TEST BY BUSCO using tetrapoda_odb10
run_BUSCO.py -i MODELS2.CDS2.faa -o BUSCO2 -l ../tetrapoda_odb10 -m proteins -c 24 -f
mv run_BUSCO2/full_table_BUSCO2.tsv ALLMODELS.BUSCO.tsv

##functional annotation using custom PROTDB
lastal -P 24 -p BL80 -m 100 -K1 PROTDB MODELS2.CDS2.faa > MODELS2.CDS2_vs_PROTDB.maf
grep ">" proteins.reference.faa | awk '{l=length($1);print substr($0,2,l-1)"\t"substr($0,l+2)}' > PROT-DB.desc.txt
maf-convert tab MODELS2.CDS2_vs_PROTDB.maf| grep -v '^#'| sed "s/_frame.\t/\t/g" | sed "s/EG2=//g" | sort --buffer-size=128G --parallel=16 -k7,7V -k13,13g | awk '{if($7!=o && $13<=1){print};o=$7}' | awk -v infile=PROT-DB.desc.txt 'BEGIN{while(getline l<infile){split(l,a,"\t");h[a[1]]=a[2]}} {print $7"\t"h[$2]" | "$2"\t"$1" e-val: "$13}' > MODELS2.CDS2_vs_PROTDB.description.txt

##functional annotation with EggNog
emapper.py --cpu 24 --mp_start_method forkserver --data_dir /data/software/HANNO/EGGNOGG-DBs -o out --output_dir ./ --temp_dir ./ --override -m diamond --dmnd_ignore_warnings -i MODELS2.CDS2.faa --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel  > emapper.out  2>emapper.err
mv out.emapper.annotations ALLMODELS.eggnog.description.txt

##Final Clean-UP
mv MODELS2.CDS2.bed12 ALLMODELS.bed12
mv MODELS2.CDS2.faa ALLMODELS.faa
mv MODELS2.CDS2.fa ALLMODELS.cds.fa
mv MODELS2.CDS2.mRNA.fa ALLMODELS.mRNA.fa
mv MODELS2.CDS2_vs_PROTDB.description.txt ALLMODELS.lastp.description.txt
#CLUSTER MODELS BY CDS
/data/software/HANNO/scripts/CLUSTERBYCDS.sh ALLMODELS.bed12 0.1
ls | grep -Ev 'asm.fa|ALLMODELS.BUSCO.tsv|ALLMODELS.bed12|ALLMODELS.faa|ALLMODELS.mRNA.fa|ALLMODELS.cds.fa|ALLMODELS.lastp.description.txt|ALLMODELS.eggnog.description.txt' | xargs rm -rf

##create Final DB that contains all results
bash /data/software/HANNO/scripts/CREATE-DB.sh ALLMODELS.bed12 ALLMODELS.bed12.model_clusters.tsv ALLMODELS.lastp.description.txt ALLMODELS.eggnog.description.txt ALLMODELS.BUSCO.tsv > ALLMODELS-FINAL.bedDB
rm -f ALLMODELS.bed12 ALLMODELS.bed12.model_clusters.tsv ALLMODELS.cds.fa ALLMODELS.eggnog.description.txt ALLMODELS.faa ALLMODELS.lastp.description.txt ALLMODELS.mRNA.fa ALLMODELS.BUSCO.tsv orf.length
sort -k14,14n -k5,5rn ALLMODELS-FINAL.bedDB | awk '{if(o!=$14){print};o=$14}' > BESTMODELS-FINAL.bedDB

OUTPUT=$(ls ../tetrapoda_odb10/hmms/ | wc -l)
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\-' | grep  -E 'Complete|Duplicated|busco1' | cut -f 1 | sort | uniq    | awk -v max=$OUTPUT 'BEGIN{i=0} {i++} END{print "\nAnalysis of BUSCOs ( n="max" ) in BESTMODELS-FINAL.bedDB\nC: n="i-1" / "100*(i-1)/max" percent"}'
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\-' | grep  -E 'Complete|Duplicated|busco1' | cut -f 1 | sort | uniq -u | awk -v max=$OUTPUT 'BEGIN{i=0} {i++} END{print "S: n="i-1" / "100*(i-1)/max" percent"}'
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\-' | grep  -E 'Complete|Duplicated|busco1' | cut -f 1 | sort | uniq -d | awk -v max=$OUTPUT 'BEGIN{i=0} {i++} END{print "D: n="i" / "100*i/max" percent"}'
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\-' | grep -Ev 'Complete|Duplicated'        | cut -f 1 | sort | uniq    | awk -v max=$OUTPUT 'BEGIN{i=0} {i++} END{print "F: n="i-1" / "100*(i-1)/max" percent"}'
cut -f 39-42 BESTMODELS-FINAL.bedDB | grep -v '\-' | grep  -E 'Complete|Duplicated|busco1|Fragmented' | cut -f 1 | sort | uniq | awk -v max=$OUTPUT 'BEGIN{i=0} {i++} END{print "M: n="max-(i-1)" / "100-100*(i-1)/max" percent\n"}'

##OUTPUT fasta files of mRNA, CDS and aminoacid
cut -f 1-12 BESTMODELS-FINAL.bedDB | bedtools getfasta -split -s -name -bed - -fi asm.fa -fo BESTMODELS-FINAL.mRNA.fa
cut -f 1-12 BESTMODELS-FINAL.bedDB | awk -f /data/software/HANNO/scripts/bed12ToGTF.awk | /data/software/HANNO/scripts/CDS_gtfToBed12 | bedtools getfasta -split -s -name -bed - -fi asm.fa -fo BESTMODELS-FINAL.CDS.fa
/data/software/HANNO/scripts/TRANSLATE.sh BESTMODELS-FINAL.CDS.fa > BESTMODELS-FINAL.AA.faa
##OUTPUT gtf
cut -f 1-12 BESTMODELS-FINAL.bedDB | awk -f /data/software/HANNO/scripts/bed12ToGTF.awk > BESTMODELS-FINAL.gtf
rm asm.fa asm.fa.fai

date
##END HANNO pipeline
