#!/bin/bash

#
# runs DESeq
#

#PBS -l walltime=#walltimeHours:00:00
#PBS -l ncpus=1
#PBS -l mem=50gb

#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe

module load R/#rVersion
module load bedtools/#bedtoolsVersion

RESULT_DIR=#resultsDir
R_SCRIPT=#rScript 
GENE_ANN=#geneAnn
PROMOTER_ANN=#promoterAnn
WINDOW=#window
DEPLOYMENT_SERVER=#deploymentServer
SUMMARY_DEPLOYMENT=#summaryDeployment
ANALYSYS=`basename $R_SCRIPT`

R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

#sort gain/loss
sed 1d $RESULT_DIR/nbinom.sig.tsv | sed 's/\"//g' | perl -e 'while(<>){@data = split(/\t/, $_); $data[1] =~ /(\S+):(\d+)-(\d+)/; print "$1\t$2\t$3\n" if $data[5] < 1 }' > $RESULT_DIR/nbinom.sig.loss.bed
sed 1d $RESULT_DIR/nbinom.sig.tsv | sed 's/\"//g' | perl -e 'while(<>){@data = split(/\t/, $_); $data[1] =~ /(\S+):(\d+)-(\d+)/; print "$1\t$2\t$3\n" if $data[5] > 1 }' > $RESULT_DIR/nbinom.sig.gain.bed

#merge
WINDOW_GAP=$(( $WINDOW + 2 ))
mergeBed -d $WINDOW_GAP -i $RESULT_DIR/nbinom.sig.loss.bed > $RESULT_DIR/nbinom.sig.loss.merge.bed
mergeBed -d $WINDOW_GAP -i $RESULT_DIR/nbinom.sig.gain.bed > $RESULT_DIR/nbinom.sig.gain.merge.bed

#return overlapping genes and promoters
windowBed -w 100 -a $RESULT_DIR/nbinom.sig.loss.bed -b $GENE_ANN |cut -f 7|sort|uniq|cut -f1 -d ';' > $RESULT_DIR/nbinom.sig.loss.genes.ID
windowBed -w 100 -a $RESULT_DIR/nbinom.sig.gain.bed -b $GENE_ANN |cut -f 7|sort|uniq|cut -f1 -d ';' > $RESULT_DIR/nbinom.sig.gain.genes.ID

windowBed -a $RESULT_DIR/nbinom.sig.loss.bed -b $PROMOTER_ANN |cut -f 7|sort|uniq|cut -f1 -d ';' > $RESULT_DIR/nbinom.sig.loss.promoters.ID
windowBed -a $RESULT_DIR/nbinom.sig.gain.bed -b $PROMOTER_ANN |cut -f 7|sort|uniq|cut -f1 -d ';' > $RESULT_DIR/nbinom.sig.gain.promoters.ID

#make an output file: location, normalized counts per sample, mean, fold change, adjusted p-value, gene name
SAMPLES=`head -n 1 $RESULT_DIR/counts.norm.tsv`

ALL_PREF="$RESULT_DIR/nbinom.sig.loss.genes $RESULT_DIR/nbinom.sig.gain.genes $RESULT_DIR/nbinom.sig.loss.promoters $RESULT_DIR/nbinom.sig.gain.promoters"
for PREF in $ALL_PREF; do

	printf "coord\t$SAMPLES\tMeanA\tMeanB\tfoldChange\tpadj\tgeneInfo\n" > $PREF.info

	EXT="${PREF##*.}"
	BED="${PREF%.*}"
	if [ $EXT == "genes" ]; then
		windowBed -w 100 -a $BED.bed -b $GENE_ANN > $PREF.bed
	else 
		windowBed -a $BED.bed -b $PROMOTER_ANN > $PREF.bed
	fi

	while read LINE; do
		LOC=`echo "$LINE"|awk -F $'\t' '{print $1":"$2"-"$3}'`
		COUNTS=`grep $LOC $RESULT_DIR/counts.norm.tsv | sed 's/\"//g'`
		INFO=`grep $LOC $RESULT_DIR/nbinom.sig.tsv | awk -F $'\t' '{print $4"\t"$5"\t"$6"\t"$9}'`
		GENE=`echo "$LINE" | awk -F $'\t' '{print $7}'`
		printf "$COUNTS\t$INFO\t$GENE\n" >> $PREF.info
	done < $PREF.bed

	if [ $EXT == "promoters" ]; then
		cat $PREF.info|uniq > $PREF.tmp
		mv $PREF.tmp $PREF.info
	fi

	rm $PREF.bed
done

cat $GENE_ANN | cut -f 4 | cut -d ';' -f 1 > $RESULT_DIR/Ensembl.gene.ID

printf "<HTML>" > $RESULT_DIR/index.html
printf "<HEAD><TITLE>DESeq - differential methylation</TITLE></HEAD>" >> $RESULT_DIR/index.html
printf "<BODY>" >> $RESULT_DIR/index.html
printf "<P><TABLE CELLPADDING=5>" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = counts.norm.tsv>counts.norm.tsv</A><TD>table of normalized counts for each sample for each ${WINDOW} bp bin" >> $RESULT_DIR/index.html

printf "<TR><TD><A HREF = dispersion_plot.png>dispersion_plot.png<TD>empirical (black dots) and fitted (red line) dispersion values plotted against the mean of the normalized counts" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = FC_plot.png>FC_plot.png<TD>scatterplot of normalized mean versus log2 fold change (red dots for differentially methylated bins with adjusted p-value below 0.1)" >> $RESULT_DIR/index.html

if [ $ANALYSYS == "DESeq.R" ]
then
	printf "<TR><TD><A HREF = PCA.png>PCA.png<TD>PCA plot for 1000 most highly methylated bins" >> $RESULT_DIR/index.html
	printf "<TR><TD><A HREF = vst_heatmap.png>vst_heatmap.png<TD>heatmap of count table for 1000 most highly methylated bins" >> $RESULT_DIR/index.html
	printf "<TR><TD><A HREF = distances_heatmap.png>distances_heatmap.png<TD>heatmap of the sample-to-sample distances for all bins with non-zero counts" >> $RESULT_DIR/index.html
	printf "<TR><TD><A HREF = filtering_param_plot.png>filtering_param_plot.png<TD>this plot helps you to choose the combination of filtering statistics and cut-off score that produce results with the highest number of differentially methylated bins after the independent filtering step" >> $RESULT_DIR/index.html
else
	printf "<TR><TD><A HREF = original_vs_filt.tsv>original_vs_filt.tsv<TD>table that compare the number of bins differentially methylated at FDR = 10% in the original analysis and after filtering bins with low counts" >> $RESULT_DIR/index.html
	printf "<TR><TD><A HREF = rank_vs_pval_plot.png>rank_vs_pval_plot.png<TD>scatterplot of bins ranked by filtering statistics vs its log10(p-value) for differential methylation" >> $RESULT_DIR/index.html
fi
printf "<TR><TD>" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.tsv>nbinom.tsv<TD>list of all ${WINDOW} bp bins with DESeq results sorted by p-value" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.tsv>nbinom.sig.tsv<TD>list of differentially methylated bins with DESeq results sorted by p-value" >> $RESULT_DIR/index.html

printf "<TR><TD><B>gain of methylation</B>" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.gain.bed>nbinom.sig.gain.bed<TD>list of differentially methylated bins (adjusted p-value below 0.1) in BED format" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.gain.merge.bed>nbinom.sig.gain.merge.bed<TD>list of merged differentially methylated bins (adjusted p-value below 0.1) in BED format (gaps of $WINDOW bp allowed)" >> $RESULT_DIR/index.html

printf "<TR><TD><A HREF = nbinom.sig.gain.genes.info>nbinom.sig.gain.genes.info<TD>list of differentially methylated bins overlapping Ensembl genes (gene body + flanking 100 bp)" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.gain.genes.ID>nbinom.sig.gain.genes.ID<TD>list of Ensembl gene IDs that overlap with differentially methylated bins" >> $RESULT_DIR/index.html

printf "<TR><TD><A HREF = nbinom.sig.gain.promoters.info>nbinom.sig.gain.promoters.info<TD>list of differentially methylated bins overlapping promoters of Ensembl genes (1500 bp upstream and downstream of TSS)" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.gain.promoters.ID>nbinom.sig.gain.promoters.ID<TD>list of Ensembl gene IDs which promoters overlap with differentially methylated bins" >> $RESULT_DIR/index.html

printf "<TR><TD><B>loss of methylation</B>" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.loss.bed>nbinom.sig.loss.bed<TD>list of differentially methylated bins (adjusted p-value below 0.1) in BED format" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.loss.merge.bed>nbinom.sig.loss.merge.bed<TD>list of merged differentially methylated bins (adjusted p-value below 0.1) in BED format (gaps of $WINDOW bp allowed)" >> $RESULT_DIR/index.html

printf "<TR><TD><A HREF = nbinom.sig.loss.genes.info>nbinom.sig.loss.genes.info<TD>list of differentially methylated bins overlapping Ensembl genes (gene body + flanking 100 bp)" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.loss.genes.ID>nbinom.sig.loss.genes.ID<TD>list of Ensembl gene IDs that overlap with differentially methylated bins" >> $RESULT_DIR/index.html

printf "<TR><TD><A HREF = nbinom.sig.loss.promoters.info>nbinom.sig.loss.promoters.info<TD>list of differentially methylated bins overlapping promoters of Ensembl genes (1500 bp upstream and downstream of TSS)" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.loss.promoters.ID>nbinom.sig.loss.promoters.ID<TD>list of Ensembl gene IDs which promoters overlap with differentially methylated bins" >> $RESULT_DIR/index.html

printf "<TR><TD><TD>" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = Ensembl.gene.ID>Ensembl.gene.ID<TD>list of Ensembl gene IDs to use as a background in DAVID analysis" >> $RESULT_DIR/index.html

scp $RESULT_DIR/* $DEPLOYMENT_SERVER:$SUMMARY_DEPLOYMENT
ssh $DEPLOYMENT_SERVER chmod 0664 $SUMMARY_DEPLOYMENT/*
