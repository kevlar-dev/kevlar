FROM debian:stable
MAINTAINER Daniel Standage <daniel.standage@gmail.com>

WORKDIR /workdir/
COPY ./ /workdir/

ENV PACKAGES locales git wget autoconf build-essential libcurl4-openssl-dev \
             zlib1g-dev libbz2-dev libbz2-dev liblzma-dev libncurses5-dev \
             python3-dev python3-pip python3-venv python3-wheel
RUN apt-get update && apt-get install -y ${PACKAGES} && apt-get clean


RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8


RUN wget https://github.com/samtools/htslib/archive/1.8.tar.gz \
    &&  tar -xzf 1.8.tar.gz && cd htslib-1.8/ \
    && autoconf && autoheader && ./configure && make && make install \
    && cd - && rm -r htslib-1.8/ 1.8.tar.gz
RUN wget https://github.com/samtools/samtools/archive/1.8.tar.gz \
    && tar -xzf 1.8.tar.gz && cd samtools-1.8/ \
    && autoconf && autoheader && ./configure && make && make install \
    && cd - && rm -r samtools-1.8/ 1.8.tar.gz
RUN wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz \
    && tar -xzf v0.7.17.tar.gz && cd bwa-0.7.17/ \
    && make && cp bwa /usr/local/bin/bwa \
    && cd - && rm -r bwa-0.7.17/ v0.7.17.tar.gz


RUN pip3 install --upgrade pip 'setuptools>=32.0.0'
RUN pip3 install cython==0.28.3 pysam==0.14.1 networkx==2.1 pandas==0.23.1 scipy==1.1.0 matplotlib==2.2.0
RUN pip3 install git+https://github.com/dib-lab/khmer.git@6c893074ea005589c230fb7cb3712f0b258f42fc
RUN pip3 install .


RUN rm -r /workdir/*

CMD bash
