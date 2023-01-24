## BUILD DOCKER AUTOMATICLY
docker build -t delta-svm:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag delta-svm:v0.0.0 ndatth/delta-svm:v0.0.0

docker run -it -v /sigma4:/sigma4 --name delta-svm ndatth/delta-svm:v0.0.0

cd bsgenome
# forge and install BSgenome of atlantic salmon v3.1
bash install.sh 
# run tandem repeat finder to find tandem repeat regions
bash tandem_repeat_finder.sh

for chr in {1..29}
do
  ./trf_dat2bed.py --dat trf_bed/${chr}.fa.2.7.7.80.10.50.12.dat --bed trf_bed/${chr}.bed
  cat trf_bed/${chr}.bed >> trf_bed/atlantic_salmon_v3.1_trf.bed
done




docker push ndatth/delta-svm:v0.0.0
echo DONE

### test docker

docker run -it --rm -v /sigma4:/sigma4 --name delta-svm ndatth/delta-svm:v0.0.0
