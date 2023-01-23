## BUILD DOCKER AUTOMATICLY
docker build -t delta-svm:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag delta-svm:v0.0.0 ndatth/delta-svm:v0.0.0
docker push ndatth/delta-svm:v0.0.0
echo DONE


### test docker

docker run -it --rm -v /sigma4:/sigma4 --name delta-svm ndatth/delta-svm:v0.0.0
