#!/bin/bash

#
# runs DESeq
#

#PBS -l walltime=#walltimeHours:00:00
#PBS -l ncpus=1
#PBS -l mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

module load R/#rVersion

RESULT_DIR=#resultsDir
R_SCRIPT=#rScript 
export GFF_PATH=#gffPath
export GO_PATH=#goPath
DEPLOYMENT_SERVER=#deploymentServer
SUMMARY_DEPLOYMENT=#summaryDeployment
ANALYSYS=`basename $R_SCRIPT`

#run R script
R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

#retrive annotation for DE genes
perl -e 'print "gene_name\tgene_id\tgene_biotype\tfold_change\tadj_pval\tGO_term\n";' > $RESULT_DIR/nbinom.sig.ann.tsv
perl -e 'while(<>){

     chomp();
     s/"//g;
     @data = split(/\t/,$_); 
     $id = $data[1];
     next unless $id =~ /ENS/;

     $name = `grep $id $ENV{'GFF_PATH'}|head -n 1|cut -f 9`;
     $name =~ s/.*gene_name \"(.*?)\".*/$1/;
     chomp($name);

     $biotype = `grep $id $ENV{'GFF_PATH'}|head -n 1|cut -f 9`;
     $biotype =~ s/.*gene_biotype \"(.*?)\".*/$1/;
     chomp($biotype);

     $go = `grep $id $ENV{'GO_PATH'}|cut -f 2,3 |grep -vP "molecular_function|cellular_component|biological_process"|tr -d ":"|tr "\t" ":"|tr "\n" ";"`;

     printf("$name\t$id\t$biotype\t%.2f\t%.2e\t$go\n",$data[5],$data[8]);

}' $RESULT_DIR/nbinom.sig.tsv >> $RESULT_DIR/nbinom.sig.ann.tsv

#sort up- and down- regulated genes with at least25% change
awk '$6 >= 1.25' $RESULT_DIR/nbinom.sig.tsv|cut -f 2| tr -d '"'|grep ENS > $RESULT_DIR/nbinom.sig.25_up.ID
awk '$6 <= 0.8' $RESULT_DIR/nbinom.sig.tsv|cut -f 2| tr -d '"'|grep ENS > $RESULT_DIR/nbinom.sig.25_down.ID
cat $RESULT_DIR/nbinom.tsv|cut -f 2|tr -d '"'|grep ENSG > $RESULT_DIR/DAVID_bg.ID

#deploy files on eliot
printf "<HTML>" > $RESULT_DIR/index.html
printf "<HEAD><TITLE>DESeq - differential gene expression</TITLE></HEAD>" >> $RESULT_DIR/index.html
printf "<BODY>" >> $RESULT_DIR/index.html
printf "<P><TABLE CELLPADDING=5>" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = counts.norm.tsv>counts.norm.tsv</A><TD>table of normalized counts for each sample for each sample" >> $RESULT_DIR/index.html

printf "<TR><TD><A HREF = dispersion_plot.png>dispersion_plot.png<TD>empirical (black dots) and fitted (red line) dispersion values plotted against the mean of the normalized counts" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = FC_plot.png>FC_plot.png<TD>scatterplot of normalized mean versus log2 fold change (red dots for DE genes with adjusted p-value below 0.1)" >> $RESULT_DIR/index.html

if [ $ANALYSYS == "DESeq.R" ]
then
	printf "<TR><TD><A HREF = PCA.png>PCA.png<TD>PCA plot for 1000 most highly expressed genes" >> $RESULT_DIR/index.html
	printf "<TR><TD><A HREF = vst_heatmap.png>vst_heatmap.png<TD>heatmap of count table for 1000 most highly expressed genes" >> $RESULT_DIR/index.html
	printf "<TR><TD><A HREF = distances_heatmap.png>distances_heatmap.png<TD>heatmap of the sample-to-sample distances for all genes with non-zero counts" >> $RESULT_DIR/index.html
	printf "<TR><TD><A HREF = filtering_param_plot.png>filtering_param_plot.png<TD>this plot helps you to choose the combination of filtering statistics and cut-off score that produce results with the highest number of DE genes after the independent filtering step" >> $RESULT_DIR/index.html
else
	printf "<TR><TD><A HREF = original_vs_filt.tsv>original_vs_filt.tsv<TD>table that compare the number of DE genes at FDR = 10% in the original analysis and after filtering genes with low counts" >> $RESULT_DIR/index.html
	printf "<TR><TD><A HREF = rank_vs_pval_plot.png>rank_vs_pval_plot.png<TD>scatterplot of genes ranked by filtering statistics vs its log10(p-value) for differential expression" >> $RESULT_DIR/index.html
fi
printf "<TR><TD>" >> $RESULT_DIR/index.html

printf "<TR><TD><A HREF = input_gc_bias.png>input_gc_bias.png<TD>lowess regression plot of gene-level counts vs gene GC-content for groups of samples. In the absence of sequencing bias regression line should be horizontal. Genome-wise distribution of genes binned by their GC-content is shown as histogram." >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = input_length_bias.png>input_length_bias.png<TD>lowess regression plot of gene-level counts vs gene length for groups of samples. In the absence of sequencing bias regression line should follow the dotted blue line (gene counts ~ gene length ). Genome-wise distribution of genes binned by their length is shown as histogram." >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = DEanalysis_gc_bias.png>DEanalysis_gc_bias.png<TD>boxplot of p-values for DE genes binned by their GC-content" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = DEanalysis_length_bias.png>DEanalysis_length_bias.png<TD>boxplot of p-values for DE genes binned by their length" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = DEgenes_gc_bias.png>DEgenes_gc_bias.png<TD>plot number of significant DE genes binned by their GC-content. In the absence of discovery bias percentage of significant DE genes across all GC-content bins should be constant." >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = DEgenes_length_bias.png>DEgenes_length_bias.png<TD>plot number of significant DE genes binned by their length. In the absence of discovery bias percentage of significant DE genes across all gene length bins should be constant." >> $RESULT_DIR/index.html
printf "<TR><TD>" >> $RESULT_DIR/index.html

printf "<TR><TD><A HREF = nbinom.tsv>nbinom.tsv<TD>list of all genes with DESeq results sorted by p-value" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.ann.tsv>nbinom.sig.ann.tsv<TD>list of significantly DE genes with annotations" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.25_up.ID>nbinom.sig.25_up.ID<TD>list of Ensembl genes with gain of expression (FC >= 1.25)" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = nbinom.sig.25_down.ID>nbinom.sig.25_down.ID<TD>list of Ensembl genes with loss of expression (FC <= 0.8)" >> $RESULT_DIR/index.html
printf "<TR><TD><A HREF = DAVID_bg.ID>DAVID_bg.ID<TD>list of Ensembl genes to use as a background in DAVID analysis" >> $RESULT_DIR/index.html

scp $RESULT_DIR/* $DEPLOYMENT_SERVER:$SUMMARY_DEPLOYMENT
ssh $DEPLOYMENT_SERVER chmod 0664 $SUMMARY_DEPLOYMENT/*
chmod 0660 $RESULT_DIR/*
