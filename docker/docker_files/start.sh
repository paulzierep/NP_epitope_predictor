
#host
#PREDICTOR_PATH="/home/paul/GitRepos/NP_epitope_predictor/tool/django_wrapper/np_epitope_predictor"
#CONDA_PATH="/home/paul/tools/anaconda3/bin/activate"

#docker
PREDICTOR_PATH=./NP_epitope_predictor/tool/django_wrapper/np_epitope_predictor/
CONDA_PATH="/opt/conda/bin/activate"

# cd ./NP_epitope_predictor/tool/django_wrapper/np_epitope_predictor/

cd $PREDICTOR_PATH

echo $(pwd)

source $CONDA_PATH

source activate NP_predictor

python manage.py runserver 0.0.0.0:8000