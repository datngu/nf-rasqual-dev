## BUILD DOCKER AUTOMATICLY
docker build -t eigenmt:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag eigenmt:v0.0.0 ndatth/eigenmt:v0.0.0
docker push ndatth/eigenmt:v0.0.0
echo DONE


### test docker

docker run -it --rm -v /sigma4:/sigma4 --name eigenmt ndatth/eigenmt:v0.0.0
docker start eigenmt
docker attach eigenmt


######## manually

docker run -it --rm -v /sigma4:/sigma4 --name eigenmt nfcore/base:2.1