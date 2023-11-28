FROM python:3.8-alpine3.16

WORKDIR /usr/src/gene_descriptions

ADD . .

RUN apk update && apk add make automake build-base gcc g++ gfortran openblas-dev python3-dev libffi-dev openjdk11 bash
RUN pip install -r requirements.txt

RUN wget https://github.com/owlcollab/owltools/releases/download/2020-04-06/owltools
RUN mv owltools /bin
RUN chmod +x /bin/owltools

CMD ["/bin/bash", "/usr/src/gene_descriptions/wormbase/exec_all_pipelines.sh"]
