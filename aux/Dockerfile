FROM ubuntu:16.04
MAINTAINER Daniel Standage <daniel.standage@gmail.com>

RUN apt-get update && apt-get install -y git build-essential python-dev python-pip bwa
RUN pip install --upgrade pip setuptools==32.0.0
RUN pip install pysam==0.11.2.2 networkx==1.11 pandas==0.20.3
RUN pip install git+https://github.com/dib-lab/khmer.git@6a1f3ec994299ceda4367d75e6aa441ebad12909

CMD bash
