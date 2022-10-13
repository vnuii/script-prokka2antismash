#!/bin/bash

###################################################################################################
# 2022-8-30
# von
# annotation-hmmer/orthofisher-bedtools-antismash pipeline
###################################################################################################

### NOTE ###
# genome文件中 ">scaffold01"，">"和后面字符间不应该有空格
# fasta文件换行符格式应都为Unix格式而不是Windows

######### 请先阅读 README文件 #########


basepath_genomeseq=/home/von/ceshi/genomeseq
base_path=/home/von/ceshi
basepath_targetprtseq=/home/von/ceshi/targetprtseq


# make name_file
cd $basepath_genomeseq
ls *.fasta >> $basepath_genomeseq/allname_genome.temp.txt 
cd $basepath_targetprtseq
ls *.fasta >> $basepath_targetprtseq/allname_prt.temp.txt
sed 's/.fasta//g' $basepath_genomeseq/allname_genome.temp.txt >> $basepath_genomeseq/allname_genome.txt
sed 's/.fasta//g' $basepath_targetprtseq/allname_prt.temp.txt >> $basepath_targetprtseq/allname_prt.txt
rm $basepath_genomeseq/allname_genome.temp.txt
rm $basepath_targetprtseq/allname_prt.temp.txt

# annotation, 需下载autoprokka.py https://github.com/stevenjdunn/autoprokka
mkdir $base_path/prokka $base_path/prtinfo $base_path/prtinfo/prtseq $base_path/prtinfo/gff
autoprokka.py -i $basepath_genomeseq -o $base_path/prokka
cp $base_path/prokka/*/*.faa $base_path/prtinfo/prtseq
cp $base_path/prokka/*/*.gff $base_path/prtinfo/gff

# convert2bed
mkdir $base_path/bed_newbed $base_path/bed_newbed/bed $base_path/bed_newbed/newbed
for i in $(<$basepath_genomeseq/allname_genome.txt);
do
convert2bed --input=gff <$base_path/prtinfo/gff/$i.gff> $base_path/bed_newbed/bed/$i.bed
done

# bed2newbed
for i in $(<$basepath_genomeseq/allname_genome.txt);
do
awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,$10}' $base_path/bed_newbed/bed/$i.bed > $base_path/bed_newbed/newbed/$i.newbed
done

# alignment
mkdir $base_path/alignment_targetprt
for i in $(<$basepath_targetprtseq/allname_prt.txt);
do 
mafft --thread 8 --auto --maxiterate 1000 $basepath_targetprtseq/$i.fasta > $base_path/alignment_targetprt/$i.aligned;
done

# hmmbuild 
mkdir $base_path/hmm $base_path/hmm/hmmbuild
for i in $(<$basepath_targetprtseq/allname_prt.txt);
do 
hmmbuild $base_path/hmm/hmmbuild/$i.hmm $base_path/alignment_targetprt/$i.aligned
done
cat $base_path/hmm/hmmbuild/*.hmm >> $base_path/hmm/hmmbuild/all.hmm

# hmmsearch
mkdir $base_path/hmm/output
for i in $(<$basepath_genomeseq/allname_genome.txt);
do 
hmmsearch -o $base_path/hmm/output/$i.hmmout $base_path/hmm/hmmbuild/all.hmm $base_path/prtinfo/prtseq/$i.faa
done

# orthofisher
mkdir $base_path/orthofisher
cd $base_path/orthofisher
cp $base_path/prtinfo/prtseq/*.faa $base_path/orthofisher
cp $base_path/hmm/hmmbuild/all.hmm $base_path/orthofisher
ls *.hmm >> hmm_name.txt
ls *.faa >> faa_name_temp.txt
awk -v FS="\t" -v OFS="\t" '{print $1,$1}' faa_name_temp.txt >> faa_name.txt
rm faa_name_temp.txt
orthofisher -m hmm_name.txt -f faa_name.txt

# uniq orthofisher_output
mkdir $base_path/orthofisher/uniq_hmmsearch
cd $base_path/orthofisher/orthofisher_output/hmmsearch_output
for i in *.hmm;
do
awk '{print $1}' $i | sort | uniq >> $base_path/orthofisher/uniq_hmmsearch/$i.uniq
done

# get genome_CDS
mkdir $basepath_genomeseq/backup1 $base_path/CDS $base_path/CDS/taget_CDS
cp $basepath_genomeseq/*.fasta $basepath_genomeseq/backup1
for i in $(<$basepath_genomeseq/allname_genome.txt);
do
bedtools getfasta -fi $basepath_genomeseq/backup1/$i.fasta -bed $base_path/bed_newbed/newbed/$i.newbed -name >> $base_path/CDS/taget_CDS/$i.cds
done

# get 外加上下游5个基因的bed文件
mkdir $base_path/bed_newbed/addbed
for i in $(<$basepath_genomeseq/allname_genome.txt);
do
	for f in $(cat $base_path/orthofisher/uniq_hmmsearch/$i.faa_all.hmm.uniq);
	do
	grep -B 5 -A 5 $f $base_path/bed_newbed/newbed/$i.newbed | head -n 1 | awk '{printf "%s\t%s\t",$1,$2}' >> $base_path/bed_newbed/addbed/$i.addbed
	grep -B 5 -A 5 $f $base_path/bed_newbed/newbed/$i.newbed | tail -n 1 | awk '{printf "%s\t%s\n",$3,$4}' >> $base_path/bed_newbed/addbed/$i.addbed
	done
done

# get RIGHT/ERROR_bed
# 由上一步获取的全部11个基因可能不在一个scallfold spade上，该步筛选出位于同一个scallfold spade上的基因作为right_bed, 不位于同一个的作为error_bed
mkdir $base_path/bed_newbed/right_addbed
mkdir $base_path/bed_newbed/error_addbed
IFS=$'\n\n'
for i in $(<$basepath_genomeseq/allname_genome.txt);
do
	for f in $base_path/bed_newbed/addbed/$i.addbed;
	do
	cat $f | awk -v FS="\t" '$2 < $3 {printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}' >> $base_path/bed_newbed/right_addbed/$i.rightbed
	cat $f | awk -v FS="\t" '$2 > $3 {printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}' >> $base_path/bed_newbed/error_addbed/$i.errorbed
	done
done

# get add_CDS
mkdir $basepath_genomeseq/backup2 $base_path/additionalCDS
cp $basepath_genomeseq/*.fasta $basepath_genomeseq/backup2
for i in $(<$basepath_genomeseq/allname_genome.txt);
do
bedtools getfasta -fi $basepath_genomeseq/backup2/$i.fasta -bed $base_path/bed_newbed/right_addbed/$i.rightbed -name >> $base_path/additionalCDS/$i.addCDS.fasta
done

# prokka_again
mkdir $base_path/additionalCDS/use_for_prokka $base_path/additionalCDS/use_for_prokka/temp $base_path/additionalCDS/use_for_prokka/use $base_path/prokka_again
for i in $(<$basepath_genomeseq/allname_genome.txt);
do
sed 's/;/ /g' $base_path/additionalCDS/$i.addCDS.fasta >> $base_path/additionalCDS/use_for_prokka/temp/temp.$i.fasta
awk '{print $1}' $base_path/additionalCDS/use_for_prokka/temp/temp.$i.fasta >> $base_path/additionalCDS/use_for_prokka/use/$i.fasta
done
autoprokka.py -i $base_path/additionalCDS/use_for_prokka/use -o $base_path/prokka_again

# antismash
mkdir $base_path/gbk_usefor_antismash $base_path/antismash_out
cp $base_path/prokka_again/*/*.gbk $base_path/gbk_usefor_antismash/
cd $base_path/gbk_usefor_antismash/
for i in *.gbk;
do
antismash $i
done
for i in $(<$basepath_genomeseq/allname_genome.txt);
do
cp -r $i $base_path/antismash_out/
done
cd $base_path
