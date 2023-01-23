#!/bin/bash
#SBATCH --ntasks=32
#SBATCH --nodes=1                
#SBATCH --job-name=TRF   
#SBATCH --mem=64G                 
#SBATCH --partition=smallmem     
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL

# Running tandem repeat finder for Atlantic salmon genome version 3.1, relased 106
#conda activate gkmsvm-salmon

echo " Running from tandem repeat finder version"
# run trf
for chr in {1..29}
do
  trf seqs/${chr}.fa 2 7 7 80 10 50 12 -d &
done
wait
rm *html


for chr in {1..29}
do
  trf_dat2bed.py --dat ${chr}.fa.2.7.7.80.10.50.12.dat --bed trf_bed/${chr}.bed
  cat trf_bed/${chr}.bed >> trf_bed/atlantic_salmon_v3.1_trf.bed
done

echo "DONE running tandem repeat finder!"

