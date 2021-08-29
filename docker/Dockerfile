FROM ubuntu:18.04

RUN apt-get update && \
	apt-get install -y build-essential git cmake autoconf libtool pkg-config

WORKDIR /genozip
RUN git clone --depth 1 https://github.com/divonlan/genozip.git /genozip && cd /genozip && make docker

WORKDIR /data
ENV PATH /genozip:$PATH