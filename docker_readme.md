# NP_epitope_predictor

## django wrapper

source start_conda.sh #activate env
python manage.py collectstatic
python manage.py runserver

## Get the conda .yml

```
source start_conda.sh #check if the django wrapper is working
conda env export --no-builds | grep -v "prefix" > requirements.yml #(inside the docker_files folder)
```

## Build the docker image

```
#change .dockerignore if needed (exclude all that is not necessary)

sudo docker build -t np_epitope_predictor:1.0 . #use --no-cache flag for new build 
sudo docker build -t np_epitope_predictor:1.0 . --no-cache
#(inside the same folder as the Dockerfile)
# must be lower case 
#or add pseudo command at position where no cache should be used
```

## Build the docker image

```
sudo docker container run -p 8000:8000 -t np_epitope_predictor:1.0
sudo docker ps #see running docker images
sudo docker kill 44decf314e47 #kill specific image
```

## Push docker image to gitlab

```
#docker build -t registry.gitlab.com/iedb-tools/nonpeptidic-predictor .
sudo docker tag np_epitope_predictor:1.0 registry.gitlab.com/iedb-tools/nonpeptidic-predictor
sudo docker login registry.gitlab.com/iedb-tools/nonpeptidic-predictor
sudo docker push registry.gitlab.com/iedb-tools/nonpeptidic-predictor
```

## Store image locally

A copy of the image is stored on the fileserver TODO[location]

```
docker save -o NP_epitope_predictor_docker_image.tar np_epitope_predictor:1.0 #store
docker load -i NP_epitope_predictor_docker_image.tar #load
```