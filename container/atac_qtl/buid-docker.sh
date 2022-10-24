## BUILD DOCKER AUTOMATICLY
docker build -t atac-qtl:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag atac-qtl:v0.0.0 ndatth/atac-qtl:v0.0.0
docker push ndatth/atac-qtl:v0.0.0
echo DONE


### test docker

docker run -it --rm -v /sigma4:/sigma4 --name qtl ndatth/atac-qtl:v0.0.0
docker start atac
docker attach atac