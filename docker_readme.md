# NP_epitope_predictor

## Get the conda .yml

```
source start_conda.sh #check if the django wrapper is working
conda env export --no-builds | grep -v "prefix" > requirements.yml
```

## Build the docker image

```
sudo docker build -t NP_epitope_predictor:1.0 . #use --no-cache flag for new build
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
sudo docker tag b6029a0c83ce registry.gitlab.com/iedb-tools/nonpeptidic-predictor
sudo docker login registry.gitlab.com/iedb-tools/nonpeptidic-predictor
sudo docker push registry.gitlab.com/iedb-tools/nonpeptidic-predictor
```