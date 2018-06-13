FROM ubuntu:18.04
MAINTAINER Daniel Standage <daniel.standage@gmail.com>

WORKDIR /src/
COPY ./ /src/

ENV PACKAGES git autoconf build-essential \
             zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev \
             python3-dev python3-venv
RUN apt-get update && apt-get install -y ${PACKAGES} && apt-get clean

RUN pip3 install --upgrade pip setuptools>=32.0.0
RUN pip3 install pysam==0.14.1 networkx==2.1 pandas==0.23.1 scipy==1.1.0 matplotlib==2.2.0
RUN pip3 install git+https://github.com/dib-lab/khmer.git@6a1f3ec994299ceda4367d75e6aa441ebad12909
RUN pip3 install .

CMD bash
