#!/bin/bash

snptest="/APP/PATH/snptest_v2.5.2/snptest_v2.5.2"
plink="/APP/PATH/plink-1.9/plink"

cd "/LOCAL/DIR" || exit 2

n="$1"

$plink --data --gen gwas/fr0207_METAG_chr$n.gen.gz --sample gwas/fr0207_METAG_chr$n.samples --oxford-single-chr $n --out chr_$n

