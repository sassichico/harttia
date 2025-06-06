#!/bin/bash
#extract and concatenate the results excluding x% as burn-in
tail -n +60001 Harttia_31_mcmc.txt > post_run1.txt
tail -n +60001 Harttia_32_mcmc.txt > post_run2.txt
tail -n +60001 Harttia_33_mcmc.txt > post_run3.txt
cat post_run1.txt post_run2.txt post_run3.txt > Harttia_combined_300SNPs.txt

tail -n +60001 Harttia_41_mcmc.txt > post_run1.txt
tail -n +60001 Harttia_42_mcmc.txt > post_run2.txt
tail -n +60001 Harttia_43_mcmc.txt > post_run3.txt
cat post_run1.txt post_run2.txt post_run3.txt > Harttia_combined_400SNPs.txt

tail -n +60001 Harttia_51_mcmc.txt > post_run1.txt
tail -n +60001 Harttia_52_mcmc.txt > post_run2.txt
tail -n +60001 Harttia_53_mcmc.txt > post_run3.txt
cat post_run1.txt post_run2.txt post_run3.txt > Harttia_combined_500SNPs.txt

#Re-add manually the header of the files, prior renumbering step

python3 fast_renumber.py Harttia_combined_300SNPs.txt Harttia_300_combinado_FINAL_mcmc.txt
python3 fast_renumber.py Harttia_combined_400SNPs.txt Harttia_400_combinado_FINAL_mcmc.txt
python3 fast_renumber.py Harttia_combined_500SNPs.txt Harttia_500_combinado_FINAL_mcmc.txt
