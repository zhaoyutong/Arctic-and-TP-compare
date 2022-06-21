#!/bin/bash
# Author: Yutao
# E-Mail: yut@im.ac.cn
# Date:2020年 12月 22日 星期二 14:23:39 CST


helpDoc(){
cat <<EOF
    Usage: bash $0 <Output dir>
EOF
}

if [ $# -lt 1 ]
then
    helpDoc
    exit 1
fi

OUT=$1
READ_PATH=/mnt/nfs/chenyy/TP_Arctic/clean_data
PHYFLASH_DB=/mnt/nfs/software/Silva_138.1_PhyloFlash_database

if [ ! -e $OUT ]
then
    mkdir -p $OUT
fi

# -lib:output file prefix only one word include "_" or "-", and the remain string not allowed
ls $READ_PATH/*_R1.fq.gz |while read i
do 
    id=$(basename $i|awk -F_ '{print $1}')
    echo  "time phyloFlash.pl -dbhome $PHYFLASH_DB  -CPUs 5 -lib $id -almosteverything -log \
        -read1 $i -read2" ${READ_PATH}/${id}_R2.fq.gz  
done

echo "Warnning:You need cd $OUT by hand!"
