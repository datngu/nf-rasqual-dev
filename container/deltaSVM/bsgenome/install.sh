# forging Atlantic salmon genome version 3.1, relased 106


# sequences downloaded from ensembl
if [ -d "seqs" ]; then
  rm seqs -f -r
fi
mkdir seqs


# tandem repeat finder 
if [ -d "trf_res" ]; then
  rm trf_bed -f -r
fi
mkdir trf_bed



echo " Downloading from Ensembl sm version"
# Download fasta files of chr 1-29
for chr in {1..29}
do
  wget "http://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/Salmo_salar.Ssal_v3.1.dna_sm.primary_assembly.${chr}.fa.gz" -O "seqs/${chr}.fa.gz"
done


echo " Extracting files"
gunzip seqs/*.gz



# build based package
cp seed_template seed
sed -i 's?PATH_TO_SEQS?'${PWD}/seqs'?g' seed
# forge
Rscript forge.R
# build
R CMD build BSgenome.Salmo.Salar.Ensembl.106
# check
#R CMD check BSgenome.Salmo.Salar.Ensembl.106_1.0.tar.gz --no-manual
# install
R CMD INSTALL BSgenome.Salmo.Salar.Ensembl.106_1.0.tar.gz 

echo "BSgenome package is INSTALLED!"




# # build mask package

# cp seed_template_mask seed_mask
# #sed -i 's?PATH_TO_SEQS?'${PWD}/seqs'?g' seed_mask
# sed -i 's?PATH_TO_TRF?'${PWD}/trf_bed'?g' seed_mask
# # forge
# Rscript forge_mask.R
# # build
# R CMD build BSgenome.Salmo.Salar.Ensembl.106.mask
# # check
# R CMD check BSgenome.Salmo.Salar.Ensembl.106.mask_1.0.tar.gz --no-manual
# # install
# R CMD INSTALL BSgenome.Salmo.Salar.Ensembl.106.mask_1.0.tar.gz 



# rm seed
# rm BSgenome.Sapiens.Ensembl.3775 -f -r
# rm BSgenome.Sapiens.Ensembl.3775.Rcheck -f -r

# Notes:
# cd `echo 'cat(system.file(package="BSgenome"))' | R --vanilla --slave`
# cd pkgtemplates/BSgenome_datapkg/
# mkdir inst
# mkdir inst/extdata

