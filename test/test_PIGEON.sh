#!/bin/bash

# This is to test the results
E=0.49

cd /PIGEON

mkdir ./test/results

# 1. GxE variance component
# TL;DR
/s/bin/python3 pigeon.py \
--gxe-var \
--gxe-sumstats ./test/WHR_GxSex_munged.sumstats.gz \
--ref-ld-chr ./1kg_ref/LDscore/LDscore. \
--w-ld-chr ./1kg_ref/weights/weights.hm3_noMHC. \
--out ./test/results/WHR_GxSex

head ./test/results/WHR_GxSex_GxE_Var.txt


# 2. Scans for PGSxE and Covariant GxE
# TL;DR
/s/bin/python3 pigeon.py \
--scan \
--gxe-sumstats ./test/WHR_GxSex_munged.sumstats.gz \
--gwas-sumstats ./test/BMI_GWAS.sumstats.gz,./test/Anorexia_GWAS.sumstats.gz \
--ref-ld-chr ./1kg_ref/LDscore/LDscore. \
--w-ld-chr ./1kg_ref/weights/weights.hm3_noMHC. \
--e-sd 0.49 \
--out ./test/results/WHR_GxSex

head ./test/results/WHR_GxSex_Scan.txt

# Conditional PGSxE
/s/bin/python3 pigeon.py \
--cond \
--gxe-sumstats ./test/WHR_GxSex_munged.sumstats.gz \
--gwas-sumstats ./test/BMI_GWAS.sumstats.gz,./test/Anorexia_GWAS.sumstats.gz \
--ref-ld-chr ./1kg_ref/LDscore/LDscore. \
--w-ld-chr ./1kg_ref/weights/weights.hm3_noMHC. \
--e-sd 0.49 \
--out ./test/results/WHR_GxSex

head ./test/results/WHR_GxSex_Cond.txt