## BUILD DOCKER AUTOMATICLY
docker build -t eigen_mt:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag eigen_mt:v0.0.0 ndatth/eigen_mt:v0.0.0
docker push ndatth/eigen_mt:v0.0.0
echo DONE


### test docker

docker run -it --rm -v /sigma4:/sigma4 --name eigen_mt ndatth/eigen_mt:v0.0.0
docker start atac
docker attach atac


######## manually

docker run -it --rm -v /sigma4:/sigma4 --name eigen_mt nfcore/base:2.1