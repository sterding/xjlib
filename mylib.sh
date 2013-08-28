#!/bin/sh
# list of bash functions

##############################
# Usage: ORFanalysis bedfile
# to get normalized CpG score of a [-1000, +1000] flanking region of TSS or middle point of input region
##############################
ORFanalysis (){
    inputbed=$1
    cat $inputbed | /home/wangj2/bin/twoBitToFa -bed=stdin /home/wangj2/nearline/2bit/mm9.2bit $inputbed.cDNA.fa
    cat $inputbed.cDNA.fa | perl getORF.pl | awk 'BEGIN{id=""; founduorf=0; }{OFS="\t";if($1!=id) {if(id!="") {if(founduorf==0) printf("%d:%d;", -1,-1); printf("\n"); founduorf=0;} id=$1;printf("%s\t%d:%d\t", $1, $2, $3);start=$2;} if($1==id && $2<start && ($2+$3-start)<50) {founduorf=1; printf("%d:%d;", $2,$3);}}END{if(founduorf==0) printf("%d:%d;", -1,-1); printf("\n");}' > $inputbed.cDNA.fa.ORF
    cat $inputbed | awk '{split($11,a,","); printf("%s\t", $4); sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} printf("%d\t", sumLen); for(i=1;i<$10;i++) printf("%d,",($6=="+")?sL[i]:(sumLen-sL[i])); print ($10==1)?"-1":"";}' > $inputbed.EEJ
    cat $inputbed.cDNA.fa | awk 'BEGIN{first=0;}{if($0~/^>/) {if(first==1) printf("\n"); print;first=1;} else printf($1); }END{printf("\n");}' | perl -slne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; printf("%s\t",$h); $found=0;while ($s=~/ATG/gi) {$found=1;printf("%d,",pos($s)-3);} if($found==0) {printf("-1,");} printf("\n");' > $inputbed.ATG
    #paste <(sort -k1,1 $inputbed.EEJ) <(sort -k1,1 $inputbed.cDNA.fa.ORF) <(sort -k1,1 $inputbed.ATG) | cut -f1-3,5-6,8 | awk '{OFS="\t"; split($6,a,",");split($4,b,":");uATG="";for(i=1;i<=length(a);i++) if(a[i]<b[1] || a[i]==-1) uATG=a[i]","uATG; else break; $6=uATG; print;}'
    paste <(sort -k1,1 $inputbed.EEJ) <(sort -k1,1 $inputbed.cDNA.fa.ORF) <(sort -k1,1 $inputbed.ATG) | cut -f1-3,5-6,8
}

##############################
# Usage: getCpGscore bedfile [tss]
# to get normalized CpG score of a [-1000, +1000] flanking region of TSS or middle point of input region
##############################
getCpGscore (){
    bedfile=$1
    option=$2
    grep -v track $bedfile | awk -v option=$option '{OFS="\t"; start=(option=="tss")?(($6=="+")?$2:$3):int(($2+$3)/2); print $1, start, start,$4,$5,$6;}' | slopBed -b 1000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | perl -ne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; @n1 = ($s=~/CG/gi); $n1=@n1; $n2 = length($s); @n3 = (($s=~/C/gi), ($s=~/G/gi)); $n3=@n3; printf("%s\t%.3f\n", $h, 4*$n2*$n1/($n3**2+1));'
}
##############################
# Usage: scanCpGisland bedfile [tss]
# to scan CpG island for a [-1000, +1000] flanking region of TSS or middle point of input region
# where the CpG island is defined by 3 conditions: (200nt windows, Obs/Exp>0.6, GC>0.5)
##############################
scanCpGisland (){
    bedfile=$1
    option=$2
    grep -v track $bedfile | awk -v option=$option '{OFS="\t"; start=(option=="tss")?(($6=="+")?$2:$3):int(($2+$3)/2); print $1, start, start,$4,$5,$6;}' | slopBed -b 1000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | perl -ne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; my @a; for($i=0;$i<(length($s)-200);$i++) {$s0=substr $s,$i,200; @n1 = ($s0=~/CG/gi); $n1=@n1; $n2 = length($s0); @n3 = ($s0=~/C/gi); @n4=($s0=~/G/gi); $n3=@n3; $n4=@n4; $obs2exp=sprintf("%.2f",($n1*$n2)/($n3*$n4+1)); $gc=sprintf("%.2f",($n3+$n4)/$n2); if($obs2exp>0.6 && $gc>0.5) {push @a, join(",",$i,$obs2exp,$gc);}} if(($#a+1)==0) {push @a, "NA";} printf("%s\t%s\n", $h, join(";",@a));'
}
##############################
# Usage: scanCpGisland bedfile [tss]
# to scan CpG island for a [-1000, +1000] flanking region of TSS or middle point of input region
# where the CpG island is defined by 3 conditions: (200nt windows, Obs/Exp>0.6, GC>0.5)
##############################
scanTATA (){
    bedfile=$1
    option=$2
    TATAmotif=/tmp/TATA/TATA.meme
    if [ ! -f $TATAmotif ];  # if not exist, then generate the motif file
    then
        # TATA motif
        [ -d /tmp/TATA/ ] || mkdir /tmp/TATA/
        echo "
        A  [16 352   3 354 268 360 222 ]
        C  [46   0  10   0   0   3   2 ]
        G  [18   2   2   5   0  20  44 ]
        T  [309  35 374  30 121   6 121 ]
        " | sed '/^$/d;s/  \[/| /;s/ \]//;s/ +/ /g'> /tmp/TATA/TATA.cm  # copied from JASPAR
        # convert to meme format, jaspar2meme from MEME Suite is required.
        jaspar2meme -cm /tmp/TATA > $TATAmotif
    fi

    if [ "$option" == "tss" ];
    then
        # scan the [-50, 0] of TSS
        grep -v track $bedfile | awk -v option=$option '{OFS="\t"; start=($6=="+")?($2-50):$3; print $1, start, start+50,$4,$5,$6;}' | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | fimo --text --verbosity 1 --output-pthresh 0.005 $TATAmotif - | awk  '{OFS="\t"; if($5=="+") print $2,$3, $8, $6}' | sort -k1,1 -k4,4gr | groupBy -g 1 -c 2,3 -o collapse,collapse
    else
        # scan the whole region
        grep -v track $bedfile | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | fimo --text --verbosity 1 --output-pthresh 0.005 $TATAmotif - | awk  '{OFS="\t"; if($5=="+" || $5=="-") print $2,$3, $8, $6}' | sort -k1,1 -k4,4gr | groupBy -g 1 -c 2,3 -o collapse,collapse
    fi
}

##############################
# Usage: getBIP bedfile1 bedfile2
# to get bidirectional promoter (BIP) regions defined by TSS in befile1 and TSS in bedfile2 (not within themselves)
# where BIP is defined as a region between two TSS with different orientation and distance <1000nt (allowing up to 100nt overlap).
# When multiple antisense transcripts are found, only report the longest bidirectional promoter.
# Output format: chr, name for Tx in bedfile1, start TSS, strand, antisense Tx name in bedfile2, antisense Tx number, antisense Tx TSS, Tx strand, length of bidirectional promoter
##############################

getBIP ()
{
    bedfile1=$1; bedfile2=$2;
    if [ "$bedfile1" == "$bedfile2" ];
    then
        intersectBed -a <(awk '{OFS="\t"; left=($6=="+")?($2-500):($3-50); if(left<0) left=0; print $1,left,left+550,$4,$5,$6}' $bedfile1) -b <(awk '{OFS="\t"; left=($6=="+")?($2-500):($3-50); if(left<0) left=0; print $1,left,left+550,$4,$5,$6}' $bedfile2) -S -wo |  awk -v up=$FLANKINGup -v down=$FLANKINGdown '{OFS="\t"; if($6=="-" && $12=="+") {s=$2+down; e=$8+up; len=(s>e)?(s-e):(e-s); print $1,$4,s,$6, $10,e,$12, len;}}' | sort -k2,2 -k8,8rn | awk 'BEGIN{ID=""}{OFS="\t"; if($2!=ID) {print; ID=$2;}}'
    else
        intersectBed -a <(awk '{OFS="\t"; left=($6=="+")?($2-500):($3-50); if(left<0) left=0; print $1,left,left+550,$4,$5,$6}' $bedfile1) -b <(awk '{OFS="\t"; left=($6=="+")?($2-500):($3-50); if(left<0) left=0; print $1,left,left+550,$4,$5,$6}' $bedfile2) -S -wo |  awk -v up=$FLANKINGup -v down=$FLANKINGdown '{OFS="\t"; if($6=="-" && $12=="+") {s=$2+down; e=$8+up;} if($6=="+" && $12=="-") {s=$2+up; e=$8+down;} len=(s>e)?(s-e):(e-s); print $1,$4,s,$6, $10,e,$12, len;}' | sort -k2,2 -k8,8rn | awk 'BEGIN{ID=""}{OFS="\t"; if($2!=ID) {print; ID=$2;}}'
    fi
}


##############################
# Usage: mapping_junction_reads mapped2splicing_rnaseq -q $readsfile 1 junctions2intron.pi.splicing.bed
# Use for smallRNA reads mapping (probably merge with the below one in the future)
# BUGFIX: replace Score with 0 (instead of $5, which can introduce duplicate intron due to the different score for same item from cufflinks and trinity) (May-16-2013)
##############################

mapping_junction_reads2 ()
{
    outputdir=$1
    filetype=$2 # -f:fasta -q:fastq
    readsfile=$3 # fa or fq
    mismatch=$4
    gene_annotation_bed=$5
    FLANKING=200; HANGOVER=3; cpu=8

    # for debug
    #outputdir=mapped2splicing_pirna; filetype="-f"; readsfile=/home/dongx/nearline/Xin/smallRNA/Phil.SRA.wt.ox.6wk.testes.raw.uniq.reads.fa.gz; mismatch=1; gene_annotation_bed=junctions2intron.splicing.bed; READS_MAPPED2GENOME_SAM=mapped2genome_rnaseq_accepted_hits_R1.sam

    # Phred score encoding
    if [ "$filetype" == "-q" ];
    then
        phred=`getphred $readsfile`
        [ "$phred" == "Phred+33" ] && scoreoption="--phred33-quals";
        [ "$phred" == "Phred+64" ] && scoreoption="--phred64-quals";
        [ "$phred" == "Solexa+64" ] && scoreoption="--solexa-quals";
        echo "fastq format with score format: $scoreoption;"
    else
        scoreoption="";
        echo "fasta format"
    fi

    ## (-200,+200) flanking region around splicing site
    [ -d $outputdir ] || mkdir $outputdir
    cat $gene_annotation_bed | awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B=""; for(i=1;i<length(a)-1;i++) {A=A""(b[i+1]-b[i]-a[i])",";B=B""(b[i]+a[i]-(b[1]+a[1]))",";} if($10>1) print $1,$2+a[1], $3-a[length(a)-1], "intron",0,$6,$2+a[1], $3-a[length(a)-1],$9,$10-1,A,B;}' | bed12ToBed6 | sort -u | awk -v fk=$FLANKING '{OFS="\t"; flanking=($2>fk)?fk:$2; print $1,$2-flanking,$3+flanking,$1"_"($2-flanking)"_"($3-$2)"_"$6, $5,$6,$2-flanking,$3+flanking,"0,0,0",2,flanking","flanking,"0,"(flanking+$3-$2);}' | twoBitToFa -bed=stdin $GENOME/mm9/Sequence/WholeGenomeFasta/genome.2bit $outputdir/minigenome.fa

    echo "...# makeing index ..."
    oldpwd=$PWD; cd $outputdir
    bowtie-build -q minigenome.fa minigenome
    cd $oldpwd;
    echo "...# mapping to junctions..."
    echo "bowtie -v $mismatch -a --best --un unmapped.fa -p $cpu $filetype -S --sam-nohead --sam-nosq genome $readsfile 2>log | sam2bed -v bed12=T -v sCol=NH > allmap.bed"
    export BOWTIE_INDEXES=$outputdir
    # Note: only take the reads mapping to forward strand (not reverse-complementary one) since we are mapping reads to transcriptome
    bowtie -v $mismatch --norc -a --best --un $outputdir/unmapped.fa -p $cpu $filetype minigenome $readsfile 2>log | awk '{OFS="\t"; print $3,$4,$4+length($5), $1,0,$2}' > $outputdir/allmap.bed

    echo "...# extract all junction reads"
    ## reads perfectly aligned to the junction, convert the coodinates, and recover the strand info from the raw junction
    awk -v fk=$FLANKING -v ho=$HANGOVER '{OFS="\t"; split($1, a, "_"); if($2<(fk-ho) && $3>(fk+ho)) print a[1], $2+a[2], a[2]+$3+a[3], $4, $5, a[4], $2+a[2], a[2]+$3+a[3], "0", 2, (fk-$2)","($3-fk), "0,"(fk-$2+a[3]);}' $outputdir/allmap.bed > $outputdir/junctions.bed

}

##############################
# Usage: mapping_junction_reads mapped2splicing_rnaseq -q $readsfile 1 junctions2intron.pi.splicing.bed
# Used for junction analysis
##############################

mapping_junction_reads ()
{
    outputdir=$1
    filetype=$2 # -f:fasta -q:fastq
    readsfile=$3 # fa or fq
    mismatch=$4
    inputbed=$5
    READS_MAPPED2GENOME_SAM=$6  #$HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/rnaseq_reads_mapped2genome.tab
    FLANKING=200; HANGOVER=2; cpu=8

    # for debug
    #outputdir=mapped2splicing_pirna; filetype="-f"; readsfile=/home/dongx/nearline/Xin/smallRNA/Phil.SRA.wt.ox.6wk.testes.raw.uniq.reads.fa.gz; mismatch=1; inputbed=junctions2intron.splicing.bed; READS_MAPPED2GENOME_SAM=mapped2genome_rnaseq_accepted_hits_R1.sam

    # Phred score encoding
    if [ "$filetype" == "-q" ];
    then
        phred=`getphred $readsfile`
        [ "$phred" == "Phred+33" ] && scoreoption="--phred33";
        [ "$phred" == "Phred+64" ] && scoreoption="--phred64";
        [ "$phred" == "Solexa+64" ] && scoreoption="--solexa-quals";
        echo "fastq format with score format: $scoreoption;"
    else
        scoreoption="";
        echo "fasta format"
    fi

    echo "...# makeing index ..."
    ## left splicing region (-200,+200) around left splicing site
    [ -d $outputdir ] || mkdir $outputdir
    export BOWTIE2_INDEXES=$GENOME/mm9/Sequence/Bowtie2Index/
    bedToGenePred $inputbed stdout | genePredToGtf file stdin $outputdir/genes.gtf
    oldpwd=$PWD
    cd $outputdir
    echo "...#  Entering $outputdir ..."
    gtf_to_fasta genes.gtf $BOWTIE2_INDEXES/genome_offrate3.fa minigenome.fa
    bowtie2-build -q minigenome.fa minigenome
    echo "...# mapping ..."
    export BOWTIE2_INDEXES=$outputdir
    echo "bowtie2 -x minigenome -p $cpu $filetype -U $readsfile -k 100 --no-unal -S accepted_hits.sam"
    bowtie2 -x minigenome -p $cpu $filetype -U $readsfile $scoreoption -k 100 --no-unal -S accepted_hits.sam

    echo "...# extract all junction reads"
    ## reads perfectly aligned to the junction
    awk -v fk=$FLANKING -v ho=$HANGOVER -v mm=$mismatch '{if($1!~/^@/ && $6~/^[0-9]+M$/ && match($0,"NM:i:[0-"mm"]\t") && $4<=(fk-ho) && $4+length($10)>=(fk+ho)) print;}' accepted_hits.sam > accepted_hits.junction.sam
    cp accepted_hits.junction.sam accepted_hits.includinggenomicmapper.sam

    if [ "$READS_MAPPED2GENOME_SAM" != "" ];
    then
        echo "...# exclude reads mapped to genome (or splicing site?) with CIGAR line xxxM ..."
        if [ ! -f $READS_MAPPED2GENOME_SAM ];
        then
            echo "...# mapping reads to the whole genome (only report the best)..."
            export BOWTIE2_INDEXES=$GENOME/mm9/Sequence/Bowtie2Index/  #(doesnot work because bowtie2 will look for current folder first)
            echo "...# running: $HOME/bin/bowtie2-2.0.0-beta6/bowtie2 -x genome -p 8 $filetype -U $readsfile -M 1 --no-unal --no-hd --no-sq -S $READS_MAPPED2GENOME_SAM"
            bowtie2 -x genome -p $cpu $filetype -U $readsfile $scoreoption -M 1 --no-hd --no-sq -S $READS_MAPPED2GENOME_SAM --un ${READS_MAPPED2GENOME_SAM/sam/unmapped.fq}
            echo "...# extracting reads unmapped to genome"
            if [ "$filetype" == "-f" ]; then    # fasta
                echo "... # fasta format"
                awk '{if(NR%2==1) print substr($1,2)}' ${READS_MAPPED2GENOME_SAM/sam/unmapped.fq} | sed 's/\/[1-2]$//' > ${READS_MAPPED2GENOME_SAM/sam/unmapped.tab}
            elif [ "$filetype" == "-q" ]; then  # fastq
                echo "... # fastq format"
                # note: for Casava <1.8, the header has /1 or /2 at the end, which does not occur in the sam file
                awk '{if(NR%4==1) print substr($1,2)}' ${READS_MAPPED2GENOME_SAM/sam/unmapped.fq} | sed 's/\/[1-2]$//' > ${READS_MAPPED2GENOME_SAM/sam/unmapped.tab}
            fi
        fi

        echo "...# intersect to get reads only mapped to junctions ..."
        cut -f1 accepted_hits.junction.sam | sort -u | sort - ${READS_MAPPED2GENOME_SAM/sam/unmapped.tab} | uniq -d > accepted_hits.junction.only.tab
        fgrep -w -f accepted_hits.junction.only.tab accepted_hits.junction.sam > accepted_hits.nogenomicmapper.sam
        mv accepted_hits.nogenomicmapper.sam accepted_hits.junction.sam
    fi

    grep '>'  minigenome.fa | sed 's/>//g;s/ /\t/g' | sort -k1,1n > minigenome.fai
    n=`wc -l minigenome.fai | awk '{print $1}'`

    echo "# get reads count"
    cut -f2-3 accepted_hits.junction.sam | awk -v n=$n '{and($1,0x10)?M[$2]++:P[$2]++;}END{OFS="\t";for(i=0;i<n;i++) print i,P[i]?P[i]:0, M[i]?M[i]:0;}' > reads.tab
    echo "# get species count"
    cut -f2-4,6 accepted_hits.junction.sam | awk '{OFS="\t";$1=and($1,0x10)?"-":"+"; print;}' | sort -u | awk -v n=$n '{($1=="-")?M[$2]++:P[$2]++;}END{OFS="\t";for(i=0;i<n;i++) print i, P[i]?P[i]:0, M[i]?M[i]:0;}' > species.tab
    paste minigenome.fai reads.tab species.tab | cut -f3,5,6,8,9 | sort -k1,1 > reads.species.tab
    echo "# Entering $oldpwd"
    cd $oldpwd
    echo "# done"
}

die ()
{
    echo $1; exit;
}
