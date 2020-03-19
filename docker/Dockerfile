FROM ubuntu:18.04

RUN apt-get update && \
	apt-get install -y build-essential git cmake autoconf libtool pkg-config

WORKDIR /genozip
RUN git clone https://github.com/divonlan/genozip.git /genozip && cd /genozip && make

WORKDIR /data
ENV PATH /genozip:$PATH