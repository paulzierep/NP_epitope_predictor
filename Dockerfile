FROM continuumio/anaconda3

WORKDIR /project

COPY docker_files/requirements.yml /project/

COPY docker_files/start.sh /project/

COPY /tool/ /project/tool/

RUN conda env create -f requirements.yml

EXPOSE 8000

CMD /bin/bash -c ./start.sh
