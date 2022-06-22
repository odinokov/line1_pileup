#! /bin/bash
# helper file to download ref genome and LINE1 coordinates

bedtools=$(which bedtools)
samtools=$(which samtools)
aria2c=$(which aria2c)

DATA='./data'
LINE1='http://l1base.charite.de/BED/hsflnil1_8438_rm.bed'
HG38='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'

function downloadFile {
    URL=$1
    filename=$(basename ${URL})
    
    if [ ! -f ${filename} ] 
    then 
        if command -v $aria2c &> /dev/null
        then
        $aria2c ${URL} -o ${filename}.tmp
        else
        wget ${URL} -O ${filename}.tmp
        fi
    mv ${filename}{.tmp,}
    fi
}

mkdir -p $DATA && cd $DATA

# get files if not exist

if [ ! -f regulatory_potent_LINE1_gh38.bed ]
then downloadFile ${LINE1} && cat $(basename ${LINE1}) \
| awk -v OFS="\t" '{if (NR!=1) print $1, $2, $3, $4, $5, $6}' \
> regulatory_potent_LINE1_hg38.bed.tmp && mv regulatory_potent_LINE1_hg38.bed{.tmp,} 
fi

if [ ! -f hg38.fa ] 
then downloadFile ${HG38} && zcat $(basename ${HG38}) > hg38.fa.tmp && mv hg38.fa{.tmp,} 
fi

# main part

$samtools faidx hg38.fa

# NB! the output sequence will be reverse complemented
$bedtools getfasta -bedOut -s -name -fi hg38.fa -bed regulatory_potent_LINE1_hg38.bed \
> regulatory_potent_LINE1_hg38.seq.bed

cd ..
