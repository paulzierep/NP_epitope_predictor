FROM continuumio/anaconda3

WORKDIR /project

COPY docker_files/requirements.yml /project/

COPY docker_files/start.sh /project/

RUN conda env create -f requirements.yml

RUN pwd

COPY /tool/ /project/tool/

EXPOSE 8000

CMD /bin/bash -c ./start.sh
