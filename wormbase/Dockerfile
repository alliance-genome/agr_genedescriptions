FROM python:3.8-slim

WORKDIR /usr/src/gene_descriptions

ADD . .

RUN apt-get update && apt-get install -y tzdata wget
ENV TZ="America/Los_Angeles"
RUN pip install -r requirements.txt

RUN wget https://github.com/owlcollab/owltools/releases/download/2020-04-06/owltools
RUN mv owltools /bin
RUN chmod +x /bin/owltools

CMD ["/bin/bash", "/usr/src/gene_descriptions/wormbase/exec_all_pipelines.sh"]
