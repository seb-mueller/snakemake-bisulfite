#makefile
#make test_1.trimmed.fq
#make test_trimmed.k$(k)stringentAS.bw
#make -n test.trimmed.k100stringentASonSly250.bw
#file=bsseq_C1_wk10_pooled_Index1_1.fq.gz
#make -f makefilebsseq -n ${file%_1.fq.gz}_1.fq
#make -f makefilebsseq -n ${file%_1.fq.gz}_1.RmDup.fq_bismark_bt2_pe.sam
#make -f makefilebsseq -n ${file%_1.fq.gz}_1.fq.gz_bismark_bt2_pe.sam
.SECONDARY:
bedGraphToBigWig=/applications/UCSC-tools/bedGraphToBigWig
bowtiepath=/applications/bowtie2/bowtie2-2.2.6/
sampath=/applications/samtools/samtools-1.3/
btpath=/applications/bedtools/bedtools-2.25.0/bin/
deduplicate=python /scripts/deduplicate.py
bismark_dir=/applications/bismark/bismark_v0.14.5/
#protocoll specific
adapters=/applications/trimmomatic/Trimmomatic-0.32/adapters/bisul_adapters.fa
#species specific
#index=/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50
index_bismarck=/home/syngenta_ftp/files/nrgene/merged_v2/
#index_C_org=/home/syngenta_ftp/files/nrgene/C_Tomato/S_habrochaites_C_pseudo_chromosomes_plus_organelles
#index_I_org=/home/syngenta_ftp/files/nrgene/I_Tomato/S_lycopersicum_I_pseudo_chromosomes_plus_organelles
#index_IC=/home/syngenta_ftp/files/nrgene/merged/genomes_I_C_organelles_merged
chromsizes=/home/syngenta_ftp/files/nrgene/merged_v2/genomes_I_C_v2_organelles_fused.chrsizes
threads=4
phred=--phred64-quals

%_1.fq : %_1.fq.gz
	gunzip -c $*_1.fq.gz > $*_1.fq
	gunzip -c $*_2.fq.gz > $*_2.fq

%_1.RmDup.fq.gz : %_1.fq
	$(deduplicate) $*_1.fq $*_2.fq
	gzip $*_1.RmDup.fq $*_2.RmDup.fq
	rm $*_1.fq $*_2.fq
	
%_fusedICv2_pe.bam: %_1.RmDup.fq.gz
	#also creates $*_fusedICv2_unmapped_reads.bam
	echo $@ >> $*_log.txt

	$(bismark_dir)/bismark $(phred) -B $*_fusedICv2 -N 1 -X 1500 --score_min L,-0.2,-0.2 --samtools_path $(sampath) --unmapped -p $(threads) --bowtie2 --path_to_bowtie $(bowtiepath) $(index_bismarck) -1 $*_1.RmDup.fq.gz -2 $*_2.RmDup.fq.gz
	$(bismark_dir)/bismark $(phred) -B $*_fusedICv2_unmapped_reads --non_directional -N 1 --score_min L,-0.2,-0.2 --samtools_path $(sampath) -p $(threads) --bowtie2 --path_to_bowtie $(bowtiepath) $(index_bismarck) $*_fusedICv2_unmapped_reads_1.fq.gz,$*_fusedICv2_unmapped_reads_2.fq.gz
	
#CHG_OB_Plantd_SLy_BSseq140808_I295_FCC4UC8ACXX_L6_Index2_1.fq.gz_bismark_bt2_pe.txt	
# %_bismark_bt2_pe.bam: %_1.fq.gz	
#%_fusedICv2_CHH.bw: %_fusedICv2_pe.bam
%_fusedICv2_all.cov2cyt: %_fusedICv2_pe.bam
	#not doing --bedgraph option since the filenames are messed up, doing this manualy
	#p stands for paired end
	${bismark_dir}/bismark_methylation_extractor -p --no_overlap --ample_memory --ignore_r2 2 --report --multicore $(threads) $?
	${bismark_dir}/bismark_methylation_extractor -s --ample_memory --report --multicore $(threads) $*_fusedICv2_unmapped_reads.bam
	
	$(eval files_CHG = CHG_OB_$*_fusedICv2_pe.txt CHG_OT_$*_fusedICv2_pe.txt CHG_OB_$*_fusedICv2_unmapped_reads.txt CHG_OT_$*_fusedICv2_unmapped_reads.txt)
	$(eval files_CHH = CHH_OB_$*_fusedICv2_pe.txt CHH_OT_$*_fusedICv2_pe.txt CHH_OB_$*_fusedICv2_unmapped_reads.txt CHH_OT_$*_fusedICv2_unmapped_reads.txt)
	$(eval files_CpG = CpG_OB_$*_fusedICv2_pe.txt CpG_OT_$*_fusedICv2_pe.txt CpG_OB_$*_fusedICv2_unmapped_reads.txt CpG_OT_$*_fusedICv2_unmapped_reads.txt)
	${bismark_dir}bismark2bedGraph --CX -o $*_fusedICv2_CHG.bedgraph $(files_CHG)
	${bismark_dir}bismark2bedGraph --CX -o $*_fusedICv2_CHH.bedgraph $(files_CHH)
	${bismark_dir}bismark2bedGraph --CX -o $*_fusedICv2_CpG.bedgraph $(files_CpG)
	${bismark_dir}bismark2bedGraph --CX -o $*_fusedICv2_all.bedgraph $(files_CpG) $(files_CHH) $(files_CHG)
	
	rm CHG_OB_$*_fusedICv2_unmapped_reads.txt
	rm CHG_OT_$*_fusedICv2_unmapped_reads.txt
	rm CpG_OB_$*_fusedICv2_unmapped_reads.txt
	rm CpG_OT_$*_fusedICv2_unmapped_reads.txt
	rm CHH_OB_$*_fusedICv2_unmapped_reads.txt
	rm CHH_OT_$*_fusedICv2_unmapped_reads.txt
	rm CHG_OB_$*_fusedICv2_pe.txt
	rm CHG_OT_$*_fusedICv2_pe.txt
	rm CpG_OB_$*_fusedICv2_pe.txt
	rm CpG_OT_$*_fusedICv2_pe.txt
	rm CHH_OB_$*_fusedICv2_pe.txt
	rm CHH_OT_$*_fusedICv2_pe.txt
	rm $*_fusedICv2_CHG.bedgraph.gz.bismark.cov.gz
	rm $*_fusedICv2_CHH.bedgraph.gz.bismark.cov.gz
	rm $*_fusedICv2_CpG.bedgraph.gz.bismark.cov.gz
	
	${bismark_dir}coverage2cytosine --genome_folder ${index_bismarck} --CX -o $*_fusedICv2_all.cov2cyt $*_fusedICv2_all.bedgraph.gz.bismark.cov.gz
	
	gunzip $*_fusedICv2_CHG.bedgraph.gz $*_fusedICv2_CpG.bedgraph.gz $*_fusedICv2_CHH.bedgraph.gz
	${bedGraphToBigWig} $*_fusedICv2_CpG.bedgraph  ${chromsizes} $*_fusedICv2_CpG.bw
	${bedGraphToBigWig} $*_fusedICv2_CHH.bedgraph  ${chromsizes} $*_fusedICv2_CHH.bw
	${bedGraphToBigWig} $*_fusedICv2_CHG.bedgraph  ${chromsizes} $*_fusedICv2_CHG.bw
	rm $*_fusedICv2_CHG.bedgraph $*_fusedICv2_CpG.bedgraph $*_fusedICv2_CHH.bedgraph
	
	#/home/sm934/bin/bam2bw.sh $? ${chromsizes} sort

#%_mergedIC_CHH_OB.bedgraph: %_mergedIC_CHH.bedgraph
	#needs ref genome in name!
# 	
	#also creates .bismark.cov.gz files
	#separating bedgraphs to 
	#${bismark_dir}bismark2bedGraph --CX $*_CHG_OB -o $*_CHG_OB.bedgraph
	#${bismark_dir}bismark2bedGraph --CX $*_CHG_OT -o $*_CHG_OT.bedgraph
	#${bismark_dir}bismark2bedGraph --CX $*_CpG_OB -o $*_CpG_OB.bedgraph
	#${bismark_dir}bismark2bedGraph --CX $*_CpG_OT -o $*_CpG_OT.bedgraph
	#${bismark_dir}bismark2bedGraph --CX $*_CHH_OB -o $*_CHH_OB.bedgraph
	#${bismark_dir}bismark2bedGraph --CX $*_CHH_OT -o $*_CHH_OT.bedgraph

	
	
	#tail -n +2 $*_CHG_OB.bedgraph | perl -lan -e 'print "$$F[0]\t$$F[1]\t$$F[2]\t-$$F[3]"' > $*_both.bedgraph

	#tail -n +2 $*_CHG_OT.bedgraph >> $*_both.bedgraph
	#cat $*_both.bedgraph | sort -k1,1 -k2,2n > $*_CHG.bedgraph
#	${bedGraphToBigWig} $*_CHG.bedgraph  ${chromsizes} $*_mergedIC_CHG.bw
#
#	tail -n +2 $*_CpG_OB.bedgraph | perl -lan -e 'print "$$F[0]\t$$F[1]\t$$F[2]\t-$$F[3]"' > $*_both.bedgraph
#	tail -n +2 $*_CpG_OT.bedgraph >> $*_both.bedgraph
#	cat $*_both.bedgraph | sort -k1,1 -k2,2n > $*_CpG.bedgraph
#	${bedGraphToBigWig} $*_CpG.bedgraph  ${chromsizes} $*_mergedIC_CpG.bw
#
#	tail -n +2 $*_CHH_OB.bedgraph | perl -lan -e 'print "$$F[0]\t$$F[1]\t$$F[2]\t-$$F[3]"' > $*_both.bedgraph
#	tail -n +2 $*_CHH_OT.bedgraph >> $*_both.bedgraph
#	cat $*_both.bedgraph | sort -k1,1 -k2,2n > $*_CHH.bedgraph
#	${bedGraphToBigWig} $*_CHH.bedgraph  ${chromsizes} $*_mergedIC_CHH.bw
#
# 	${bedGraphToBigWig} CHG_OT_Plantd_SLy_BSseq140808_I295_FCC4UC8ACXX_L6_Index2_1.fq.gz_bismark_bt2_pe.bedgraph  ${chromsizes} CHG_OT_Plantd_SLy_BSseq140808_I295_FCC4UC8ACXX_L6_Index2_1.fq.gz_bismark_bt2_pe.bw
# 	
# 	${bedGraphToBigWig} CHG_OB_Plantd_SLy_BSseq140808_I295_FCC4UC8ACXX_L6_Index2_1.fq.gz_bismark_bt2_pe.bedgraph  ${chromsizes} CHG_OB_Plantd_SLy_BSseq140808_I295_FCC4UC8ACXX_L6_Index2_1.fq.gz_bismark_bt2_pe.bw
# 	
# 	
# 	${sampath}/samtools view -bS Plantd_SLy_BSseq140808_I295_FCC4UC8ACXX_L6_Index2_1.fq.gz_bismark_bt2_pe.sam > Plantd_SLy_BSseq140808_I295_FCC4UC8ACXX_L6_Index2_1.fq.gz_bismark_bt2_pe.bam
# 	

# %.trimmed.k$(k)stringentASonSly250.bedgraph : %.trimmed.k$(k)stringentASonSly250.bam
# 		samtools view -h $< | awk -F"\t" ' { if (/^@/ && substr($$2, 3, 1)==":") {print} else if ($$2 == 99 || $$2 == 147) {print}} ' | samtools view -bS - | $(btpath)/genomeCoverageBed -bg -ibam stdin -g $(chromsizes) > plus.bedgraph
# 		
# 		samtools view -h $< | awk -F"\t" ' { if (/^@/ && substr($$2, 3, 1)==":") {print} else if ($$2 == 83 || $$2 == 163) {print}} ' | samtools view -bS - | $(btpath)/genomeCoverageBed -bg -ibam stdin -g $(chromsizes) | perl -lane 'print "$$F[0]\t$$F[1]\t$$F[2]\t-$$F[3]\t"' > minus.bedgraph
# 		
# 		cat plus.bedgraph minus.bedgraph | sort -k1,1 -k2,2n > $@
# 		
# 		#$(btpath)/genomeCoverageBed -bg -ibam $< -g $(chromsizes) > $*.trimmed.k$(k)stringentASonSly250.minus.bedgraph
# 
# %.trimmed.k$(k)stringentASonSly250.plus.bw : %.trimmed.k$(k)stringentASonSly250.plus.bedgraph
# 	$(bedGraphToBigWig) $< $(chromsizes) $@
# 	$(bedGraphToBigWig) $*.trimmed.k$(k)stringentASonSly250.minus.bedgraph $(chromsizes) $*.trimmed.k$(k)stringentASonSly250.minus.bw
# 	#-rm $<
# 	#-rm *.trimmed.k$(k)stringentASonSly250.minus.bedgraph
# 	#-rm $*_*.fq
# 	-rm $**.sam
	
# PATH=$bowtiepath:$sampath:$PATH
# $(tophatpath)tophat -p 16 -o $base --mate-inner-dist 20 --mate-std-dev 100 --read-mismatches 2 --min-intron-length 40 --max-intron-length 10000 --max-multihits 500  $index ${base}_1.trimmed.fq.gz  ${base}_2.trimmed.fq.gz,${base}_1.trimmed.unpaired.fq.gz,${base}_2.trimmed.unpaired.fq.gz

# $(btpath)/genomeCoverageBed -bg -ibam $<  -strand + -g $(chromsizes) > $@
# $(btpath)/genomeCoverageBed -bg -ibam $<  -strand - -g $(chromsizes) > $*.trimmed.k$(k)stringentASonSly250.minus.bedgraph
#

#$(bowtiepath)/bowtie2 -q -k $(k) --score-min L,-0.6,-0.4 -x $(index) -1 $*_1.trimmed.fq.gz -2 $*_2.trimmed.fq.gz -U $*_1.trimmed.unpaired.fq.gz,$*_2.trimmed.unpaired.fq.gz -S $@ -p $(threads) >> log.txt 2>&1
# %.trimmed.k$(k)stringentASonSly250.bam : %.trimmed.k$(k)stringentASonSly250.sam
# 	$(sampath)/samtools view -bS $< > temp_unsorted.bam
#	rm $base"_unsorted".sam
#	$(sampath)/samtools sort temp_unsorted.bam $*.trimmed.k$(k)stringentASonSly250
#	rm temp_unsorted.bam
#	$(sampath)/samtools index $@


#PATH=$bowtiepath:$PATH
#/applications/tophat/tophat-2.0.12.Linux_x86_64/tophat -p 12 -o $base --max-multihits 500  $index $base".fastq.gz" 
##
# index=/data/public_data/tomato/assembly_build_2.50/S_lycopersicum_chromosomes.2.50
# bowtiepath=/applications/bowtie2/bowtie2-2.2.3/
# sampath=/applications/samtools-0.1.19/
# btpath=/applications/bedtools-2.17.0/bin/
# chromsizes=/data/public_data/tomato/assembly_build_2.50/tomato250.chrom.sizes.txt
# threads=12
# k=100


# for file in $(ls {H3,INPUT}*Sly_*_1.fq.gz); do
#file="H3_Lane1_BioRep1_TechRep1_M82_SL2239_24-Apr_Cycles16_140529_I338_FCC4M63ACXX_L6_index19_1.fq.gz"
# 	base=${file%_1.fq.gz}
# 	pair1=$base"_1"
# 	pair2=$base"_2"
# 	if [ ! -e mapping/$pair1.fq ]; then cp $pair1.fq.gz $pair2.fq.gz mapping; fi
# 	cd mapping
# 	if [ ! -e $pair1.fq ]; then gunzip $pair1.fq.gz $pair2.fq.gz; fi
# 	if [ ! -e "Rm_dupPE_"$pair1.fq.gz ]; then python /home/sm934/bin/deduplicate.py $pair1.fq $pair2.fq
# 	pair1="Rm_dupPE_"$pair1
# 	pair2="Rm_dupPE_"$pair2
# 	base=$base""$k"stringentAS"
# 	if [ ! -e $pair2.trimmed.unpaired.fq.gz ]; java -jar /applications/trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads $threads -trimlog trim.log $pair1.fq $pair2.fq $pair1.trimmed.fq $pair1.trimmed.unpaired.fq $pair2.trimmed.fq $pair2.trimmed.unpaired.fq ILLUMINACLIP:/applications/trimmomatic/Trimmomatic-0.32/adapters/scriptseq_adaptor.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> trim.stat 2>&1
# 	echo $base >> log.txt
# 	$bowtiepath/bowtie2 -q -k $k --score-min L,-0.6,-0.4 -x $index -1 $pair1.trimmed.fq -2 $pair2.trimmed.fq -U $pair1.trimmed.unpaired.fq,$pair2.trimmed.unpaired.fq -S $base"_unsorted".sam -p $threads  >> log.txt 2>&1
# 	$sampath/samtools view -bS $base"_unsorted".sam > $base"_unsorted".bam
# 	rm $base"_unsorted".sam
# 	$sampath/samtools sort $base"_unsorted".bam $base
# 	rm $base"_unsorted".bam
# 	$sampath/samtools index $base.bam
# 	#remove PCR duplicates
# 	#$sampath/samtools rmdup $base.bam $base.rmdup.bam
# 	#$sampath/samtools index $base.rmdup.bam
# 	#make uniqui
# 	$sampath/samtools view -bq 1 $base.bam > $base.unique.bam
# 	$sampath/samtools index $base.unique.bam
# 	$btpath/genomeCoverageBed -bg -ibam $base.bam -g $chromsizes > $base.bedgraph
# 	/applications/bedGraphToBigWig $base.bedgraph $chromsizes $base.bw
# 	gzip Rm*.fq 
# 	#rm ${base%}*.fq
# 	cd ..
# done
