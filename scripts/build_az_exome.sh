#!/bin/sh

#Ensembl.gtf
#1       miRNA   gene    5875931 5876019 .       -       .       gene_id "ENSG00000216045"; gene_name "AL356693.1"; gene_source "ensembl"; gene_biotype "miRNA";
#1       miRNA   transcript      5875931 5876019 .       -       .       gene_id "ENSG00000216045"; transcript_id "ENST00000401226"; gene_name "AL356693.1"; gene_source "ensembl"
#1       miRNA   exon    5875931 5876019 .       -       .       gene_id "ENSG00000216045"; transcript_id "ENST00000401226"; exon_number "1"; gene_name "AL356693.1"; gene_source
#
#subtruction
#chr1    5875992 5876019 AL356693.1      chr1    5875931 5876019 miRNA   gene    27
#chr1    5875992 5876019 AL356693.1      chr1    5875931 5876019 miRNA   transcript      27
#chr1    5875992 5876019 AL356693.1      chr1    5875931 5876019 miRNA   exon    27



grep  Ensembl.transcripts.gtf > Ensembl.transcripts.protein_coding.gtf

grep -wf gene_list.txt /ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/Ensembl.gtf | grep 'gene' | cut -f1,4,5 | bedtools slop -l 1 -r 0 -i - -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome | grch37_to_hg19.py  > genes.bed




###############
# TRANSCRIPTS #
###############

# get the union.bed
grep -v '^#' /ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/Ensembl.gtf | grep -w 'transcript' > Ensembl.transcripts.gtf
grep 'protein_coding\|nonsense_mediated_decay\|miRNA\|IG_C_gene\|IG_D_gene\|IG_J_gene\|IG_V_gene\|TR_C_gene\|TR_D_gene\|TR_J_gene\|TR_V_gene' Ensembl.transcripts.gtf > Ensembl.transcripts.protein_coding.gtf
cut -f1,4,5 Ensembl.transcripts.protein_coding.gtf > transcripts.txt
bedtools slop -l 1 -r 0 -i transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > transcripts.bed
grch37_to_hg19.py transcripts.bed | grep -v "^chrH" | grep -v "^chrG" > transcripts.fixedchr.bed
cd ..
cat padded_panels/*.bed transcripts.bed > joined.bed
sort -k1,1 -k2,2n joined.bed > joined.sorted.bed
bedtools merge -i joined.sorted.bed > union.bed

echo 'Size of the bed file'
cat union.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

# subtruct
bedtools subtract -b union.bed -a ../AZ-Exome.bed > subtraction.bed

# get the bed for annotation
awk '{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 }' Ensembl.transcripts.gtf > Ensembl.transcripts.txt
bedtools slop -l 1 -r 0 -i Ensembl.transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.transcripts.bed
grch37_to_hg19.py Ensembl.transcripts.bed | grep -v "^H" | grep -v "^G" > Ensembl.fixedchr.transcripts.bed

# intersect
bedtools intersect -a subtraction.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.annotated.bed
cut -f8 subtraction.annotated.bed | sort -u

bedtools intersect -a subtraction.back.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.back.annotated.bed
cut -f7 subtraction.back.annotated.bed | sort -u



#BIOTYPES='protein_coding\|nonsense_mediated_decay\|miRNA\|IG_C_gene\|IG_D_gene\|IG_J_gene\|IG_V_gene\|TR_C_gene\|TR_D_gene\|TR_J_gene\|TR_V_gene'
#
## get the union.bed
#grep -v '^#' /ngs/reference_data/genomes/Hsapiens/hg19/bed.20150528/Exons/Ensembl.gtf | grep -w 'transcript' > Ensembl.transcripts.gtf
#grep '${BIOTYPES}' Ensembl.transcripts.gtf > Ensembl.transcripts.protein_coding.gtf
#cut -f1,4,5 Ensembl.transcripts.protein_coding.gtf > transcripts.txt
#bedtools slop -l 1 -r 0 -i transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > transcripts.bed
#grch37_to_hg19.py transcripts.bed | grep -v "^chrH" | grep -v "^chrG" > transcripts.fixedchr.bed
#cd ..
#cat padded_panels/*.bed transcripts.bed > joined.bed
#sort -k1,1 -k2,2n joined.bed > joined.sorted.bed
#bedtools merge -i joined.sorted.bed > union.bed
#
## subtruct
#bedtools subtract -b union.bed -a ../AZ-Exome.bed > subtraction.bed
#
## get the bed for annotation
#awk '{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 }' Ensembl.transcripts.gtf > Ensembl.transcripts.txt
#bedtools slop -l 1 -r 0 -i Ensembl.transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.transcripts.bed
#grch37_to_hg19.py Ensembl.transcripts.bed | grep -v "^H" | grep -v "^G" > Ensembl.fixedchr.transcripts.bed
#
## intersect
#bedtools intersect -a subtraction.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.annotated.bed
#cut -f8 subtraction.annotated.bed | sort -u
#
#bedtools intersect -a subtraction.back.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.back.annotated.bed
#cut -f7 subtraction.back.annotated.bed | sort -u






####################
# UTR + PADDED CDS #
####################

function gtf2bed {
    cat $1 | cut -f1,4,5 | bedtools slop -l 1 -r 0 -i - -g ${GENOME} | grch37_to_hg19.py | grep -v "^chrH" | grep -v "^chrG"
}
function gtf2fullbed {
    cat $1 | bedtools slop -l 1 -r 0 -i - -g ${GENOME} | grch37_to_hg19.py | grep -v "^chrH" | grep -v "^chrG"
}
function unionbed {
    sort -k1,1 -k2,2n $1 | bedtools merge -i -
}
function bedsize {
    awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $1
}

GENOME=/ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome
ENSEMBL=/ngs/reference_data/genomes/Hsapiens/hg19/bed.20150528/Exons/Ensembl.gtf
BIOTYPES='protein_coding\|nonsense_mediated_decay\|non_stop_decay\|miRNA\|IG_C_gene\|IG_D_gene\|IG_J_gene\|IG_V_gene\|TR_C_gene\|TR_D_gene\|TR_J_gene\|TR_V_gene'
grep -v '^#' ${ENSEMBL} | grep "${BIOTYPES}" > Ensembl.exons_utr_cds.gtf
grep $'miRNA\texon' Ensembl.exons_utr_cds.gtf | gtf2bed - > Ensembl.mirna.bed
grep -w 'CDS' Ensembl.exons_utr_cds.gtf | gtf2bed - | bedtools slop -b 50 -i - -g ${GENOME} > Ensembl.csd.bed
grep -w 'UTR\|stop_codon' Ensembl.exons_utr_cds.gtf | gtf2bed - | bedtools slop -b 50 -i - -g ${GENOME} > Ensembl.utr.bed
cat Ensembl.mirna.bed Ensembl.csd.bed Ensembl.utr.bed > Ensembl.cds_utr_mirna.bed
cat ../padded_panels/*.bed ../transcripts/Ensembl.cds_utr_mirna.bed | union - > ../union.cds_utr_mirna.bed
echo 'Size of the bed file'
bedsize ../union.cds_utr_mirna.bed

##############
# SUBTRUCION #
# subtruct
bedtools subtract -b union.cds_utr_mirna.bed -a ../AZ-Exome.bed > subtraction.bed

# get the bed for annotation
awk '{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 }' Ensembl.transcripts.gtf > Ensembl.transcripts.txt
bedtools slop -l 1 -r 0 -i Ensembl.transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.transcripts.bed
grch37_to_hg19.py Ensembl.transcripts.bed | grep -v "^H" | grep -v "^G" > Ensembl.fixedchr.transcripts.bed

# intersect
bedtools intersect -a subtraction.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.annotated.bed
cut -f8 subtraction.annotated.bed | sort -u

bedtools intersect -a subtraction.back.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.back.annotated.bed
cut -f7 subtraction.back.annotated.bed | sort -u



awk '{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 }' Ensembl.exons_utr_cds.gtf | gtf2fullbed > Ensembl.exons_utr_cds.bed
bedtools intersect -a ../refseq/gaps2.bed -b Ensembl.exons_utr_cds.bed -wao > gaps2.annotated.bed

# adding refseq
bedtools slop -b 50 -g ${GENOME} -i refseq/cds.bed > refseq/cds.padded50.bed
cat union.cds_utr_mirna.bed refseq/cds.padded50.cut3.bed | unionbed > union.cds_utr_mirna_ensembl_refseq.bed


##### BASED ON EXONS ####
#
## get the union.exons.bed
#C_BIOTYPES='protein_coding\|nonsense_mediated_decay\|miRNA\|IG_C_gene\|IG_D_gene\|IG_J_gene\|IG_V_gene\|TR_C_gene\|TR_D_gene\|TR_J_gene\|TR_V_gene'
#NC_BIOTYPES='miRNA'
#ENSEMBLE=/ngs/reference_data/genomes/Hsapiens/hg19/bed.20150528/Exons/Ensembl.gtf
#grep -v '^#' ${ENSEMBLE} | grep -w 'CDS\|stop_codon' | grep '${BIOTYPES}' > Ensembl.cds.gtf
#grep -v '^#' ${ENSEMBLE} | grep -w 'UTR' | grep '${BIOTYPES}' > Ensembl.utr.gtf
#grep -v '^#' ${ENSEMBLE} | grep -w 'exon' | grep '${BIOTYPES}' > Ensembl.exon.gtf
#grep -w 'CDS' Ensembl.cds_utr.gtf | cut -f1,4,5 > Ensembl.csd.prebed
#grep -vw 'CDS' Ensembl.cds_utr.gtf | cut -f1,4,5 > Ensembl.utr.prebed
#bedtools slop -l 1 -r 0 -i Ensembl.csd.prebed -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.csd.bed
#bedtools slop -l 1 -r 0 -i Ensembl.utr.prebed -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.utr.bed
#grch37_to_hg19.py Ensembl.csd.bed | grep -v "^chrH" | grep -v "^chrG" > Ensembl.csd.fixedchr.bed
#grch37_to_hg19.py Ensembl.utr.bed | grep -v "^chrH" | grep -v "^chrG" > Ensembl.utr.fixedchr.bed
#bedtools slop -b 50 -i Ensembl.csd.fixedchr.bed -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.csd.fixedchr.slop50.bed
#cat Ensembl.utr.bed Ensembl.csd.fixedchr.slop50.bed > Ensembl.cds_utr.bed
#cd ..
#cat padded_panels/*.bed transcripts/Ensembl.cds_utr.bed > joined.cds_utr.bed
#sort -k1,1 -k2,2n joined.cds_utr.bed > joined.cds_utr.sorted.bed
#bedtools merge -i joined.cds_utr.sorted.bed > union.cds_utr.bed
#
## subtruct
#bedtools subtract -b union.bed -a ../AZ-Exome.bed > subtraction.bed
#
## get the bed for annotation
#awk '{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 }' Ensembl.transcripts.gtf > Ensembl.transcripts.txt
#bedtools slop -l 1 -r 0 -i Ensembl.transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.transcripts.bed
#grch37_to_hg19.py Ensembl.transcripts.bed | grep -v "^H" | grep -v "^G" > Ensembl.fixedchr.transcripts.bed
#
## intersect
#bedtools intersect -a subtraction.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.annotated.bed
#cut -f8 subtraction.annotated.bed | sort -u
#
#bedtools intersect -a subtraction.back.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.back.annotated.bed
#cut -f7 subtraction.back.annotated.bed | sort -u




#### RefSeq GAPS ###
# get gaps
bedtools subtract -a refseq/cds.bed -b union.cds_utr_mirna.bed > refseq/gaps.bed

# get file for annotation
grep -v '^#' ${ENSEMBL} | grep -w 'CDS\|UTR\|stop_codon' | grep ${BIOTYPES} > Ensembl.cds_utr.gtf
awk '{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 }' Ensembl.transcripts.gtf > Ensembl.transcripts.txt
bedtools slop -l 1 -r 0 -i Ensembl.transcripts.txt -g /ngs/reference_data/genomes/Hsapiens/hg19/genome/human.hg19.genome > Ensembl.transcripts.bed
grch37_to_hg19.py Ensembl.transcripts.bed | grep -v "^H" | grep -v "^G" > Ensembl.fixedchr.transcripts.bed

# intersect
bedtools intersect -a subtraction.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.annotated.bed
cut -f8 subtraction.annotated.bed | sort -u

bedtools intersect -a subtraction.back.bed -b Ensembl.fixedchr.transcripts.bed -wao > subtraction.back.annotated.bed
cut -f7 subtraction.back.annotated.bed | sort -u














