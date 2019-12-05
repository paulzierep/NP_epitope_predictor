FROM continuumio/anaconda3

WORKDIR /project

COPY docker_files/requirements.yml /project/

COPY docker_files/start.sh /project/

RUN conda env create -f requirements.yml

RUN apt install git

RUN apt install unzip

RUN git clone https://paulzierep:7722BBBBoma@github.com/paulzierep/NP_epitope_predictor.git

RUN unzip ./NP_epitope_predictor/tool/ML_data.zip -d ./NP_epitope_predictor/tool/

EXPOSE 8000

CMD /bin/bash -c ./start.sh
