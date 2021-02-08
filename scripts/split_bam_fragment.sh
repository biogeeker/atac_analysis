#!/bin/bash
samtools view -h $1 | \
  awk -v LEN=$2 '{if ($9 <= LEN && $9 >= -(LEN) && $9 != 0 || $1 ~ /^@/) print $0}' | \
  samtools view -bh -o split_fragment.bam -

## Example to get fragments of 147 bp or smaller:
#./code.sh in.bam 147
