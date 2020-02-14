FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

ENV GEM_MAPPER_VERSION=v3.3.0
ENV INSTALL_BASE=/software/opt/gem3-mapper
ENV SRC_BASE=/software/src/gem3-mapper

RUN apt-get -y update && \
    apt-get -y install \
        gcc \
	git \
        make && \
    apt-get -y clean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt-get/lists/* 

RUN mkdir -p ${SRC_BASE} && \
    cd ${SRC_BASE} && \
    git clone --recursive https://github.com/smarco/gem3-mapper.git -b ${GEM_MAPPER_VERSION} ./ && \
    chmod +x configure && \
    ./configure && \
    make

RUN ln -s ${SRC_BASE}/bin/* /usr/local/bin/

