#!/bin/bash

mkdir trim
mkdir merge
echo -e "sample-id"'\t'"absolute-filepath" > manifest-se	

# cut primer
for i in `cut -f1 sample-metadata.tsv | grep -v "#"`;do
name1=`ls rawdata | grep R1.fq.gz | grep ${i}`
cutadapt -g CCTAYGGGNBGCWSCAG -q 20 -o ${i}.R1.trim.fq rawdata/${name1}

name2=`ls rawdata | grep R2.fq.gz | grep ${i}`
cutadapt -g GACTACNNGGGTWTCTAAT -q 20 -o ${i}.R2.trim.fq rawdata/${name2}
done

# vsearch merge
for i in `cut -f1 sample-metadata.tsv | grep -v "#"`;do
R1=`ls trim | grep R1.trim.fq | grep ${i}`
R2=`ls trim | grep R2.trim.fq | grep ${i}`
vsearch --fastq_mergepairs trim/${i}.R1.trim.fq --reverse trim/${i}.R2.trim.fq --fastqout merge/${i}.fq.gz
echo -e "${i}"'\t''$PWD/merge/'"${i}.fq.gz" >> manifest-se
done

# import data
qiime tools import  --type 'SampleData[SequencesWithQuality]'  --input-path manifest-se  --output-path demux.qza  --input-format SingleEndFastqManifestPhred33V2
qiime demux summarize  --i-data demux.qza  --o-visualization demux.qzv

# denoise
qiime dada2 denoise-single  --i-demultiplexed-seqs demux.qza  --p-trunc-len 0  --o-representative-sequences rep-seqs.qza  --o-table table.qza  --o-denoising-stats dada2-stats.qza  --p-n-threads 0
## denoise data visualization
qiime metadata tabulate  --m-input-file dada2-stats.qza  --o-visualization dada2-stats.qzv
## filter features
qiime feature-table filter-features --i-table table.qza --p-min-frequency 2 --o-filtered-table table.qza
## summarize feature-table
qiime feature-table summarize  --i-table table.qza  --o-visualization table.qzv  --m-sample-metadata-file sample-metadata.tsv


mkdir -p backup/visual backup/core-metrics backup/picrust-core-metrics

SUFFIX=`pwd`

# taxa alignment
qiime feature-classifier classify-sklearn  --i-classifier classifiergg138.qza  --p-n-jobs -1  --i-reads rep-seqs.qza  --o-classification taxonomygg138.qza
qiime taxa barplot  --i-table table.qza  --i-taxonomy taxonomygg138.qza  --m-metadata-file sample-metadata.tsv  --o-visualization backup/visual/taxa-bar-plotsgg138.qzv

for i in `seq 2 7`;do
# collapse table.qza to level
qiime taxa collapse --i-table table.qza --o-collapsed-table col-table-$i.qza --i-taxonomy taxonomygg138.qza  --p-level $i
# calculate relative abundance
qiime feature-table relative-frequency  --i-table col-table-$i.qza  --o-relative-frequency-table relab-table-$i.qza
# export biom
qiime tools export --input-path relab-table-$i.qza  --output-path col-relab-$i
# convert to tsv
biom convert  -i col-relab-$i/feature-table.biom -o backup/visual/gg138-level-$i-relab.tsv  --to-tsv
done
#backup
mv taxonomygg138.qza backup/taxonomygg138-${SUFFIX##*/}.qza
trash col-relab* col-table* relab* classifier*


# diversity analysis
## construct tree (aligned)
qiime phylogeny align-to-tree-mafft-fasttree  --i-sequences rep-seqs.qza  --o-alignment aligned-rep-seqs.qza  --o-masked-alignment masked-aligned-rep-seqs.qza  --o-tree unrooted-tree.qza  --o-rooted-tree rooted-tree.qza  --p-n-threads 8
#clean
rm aligned-rep-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza 

## alpha rarefaction
qiime diversity alpha-rarefaction  --i-table table.qza  --i-phylogeny rooted-tree.qza  --p-max-depth 17209  --m-metadata-file sample-metadata.tsv  --o-visualization backup/core-metrics/gg138-alpha-rarefaction.qzv
## diversity results
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 17209  --m-metadata-file sample-metadata.tsv  --output-dir core-metrics-results-gg138

## alpha diversity
for alphaname in {evenness,shannon,faith_pd,observed_features};do
### alpha significance
qiime diversity alpha-group-significance  --i-alpha-diversity core-metrics-results-gg138/${alphaname}_vector.qza  --m-metadata-file sample-metadata.tsv  --o-visualization backup/core-metrics/gg138-alpha-${alphaname}-.qzv
done

## beta diversity
for betaname in {bray_curtis,jaccard,unweighted_unifrac,weighted_unifrac};do
### group distance significance 
qiime diversity beta-group-significance  --i-distance-matrix core-metrics-results-gg138/${betaname}_distance_matrix.qza  --m-metadata-file sample-metadata.tsv  --m-metadata-column Group  --o-visualization core-metrics-results-gg138/gg138-${betaname}_Group_significance.qzv  --p-pairwise
done

#backup
for emperor in {bray_curtis_emperor.qzv,jaccard_emperor.qzv,unweighted_unifrac_emperor.qzv,weighted_unifrac_emperor.qzv};do
mv core-metrics-results-gg138/${emperor} backup/core-metrics/gg138-beta-${emperor}
done
mv core-metrics-results-gg138/*.qzv backup/core-metrics

# picrust with denoising table and rep-seqs 
qiime picrust2 full-pipeline  --i-table table.qza  --i-seq rep-seqs.qza  --output-dir picrust2-output  --p-threads 8  --p-hsp-method mp  --p-max-nsti 2
## visualization
qiime feature-table summarize --i-table picrust2-output/pathway_abundance.qza --o-visualization func-pathway-abundance.qzv
## export data
qiime tools export --input-path picrust2-output/pathway_abundance.qza  --output-path pathabun-exported
## convert to tsv
biom convert -i pathabun-exported/feature-table.biom  -o func-pathway-abundance-${SUFFIX##*/}.tsv  --to-tsv
## clean
trash pathabun-exported

# diversity analysis
qiime diversity core-metrics --i-table picrust2-output/pathway_abundance.qza --p-sampling-depth 1007548 --m-metadata-file sample-metadata.tsv --output-dir pathabun-core-metrics-out --p-n-jobs 8
## alpha diversity
for alphaname in {evenness,shannon,observed_otus};do
### alpha significance
qiime diversity alpha-group-significance  --i-alpha-diversity pathabun-core-metrics-out/${alphaname}_vector.qza  --m-metadata-file sample-metadata.tsv  --o-visualization pathabun-core-metrics-out/pathabun-alpha-${alphaname}-.qzv
done

## beta diversity
for betaname in {bray_curtis,jaccard};do
### group distance significance 
qiime diversity beta-group-significance  --i-distance-matrix pathabun-core-metrics-out/${betaname}_distance_matrix.qza  --m-metadata-file sample-metadata.tsv  --m-metadata-column Group  --o-visualization pathabun-core-metrics-out/pathabun-beta-${betaname}_Group_significance.qzv  --p-pairwise
done

#backup
for emperor in {bray_curtis_emperor.qzv,jaccard_emperor.qzv};do
mv pathabun-core-metrics-out/${emperor} backup/picrust-core-metrics/pathabun-beta-${emperor}
done
mv pathabun-core-metrics-out/*.qzv backup/picrust-core-metrics 