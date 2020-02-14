FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

ARG GEM_MAPPER_VERSION=master
ARG INSTALL_BASE=/software/opt/gem3-mapper
ARG SRC_BASE=/software/src/gem3-mapper

RUN apt-get -y update && \
    apt-get -y install \
        gcc \
	git \
        make && \
    apt-get -y clean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt-get/lists/* 

RUN mkdir -p ${SRC_BASE} && \
    mkdir -p ${INSTALL_BASE} && \
    cd ${SRC_BASE} && \
    git clone --recursive https://github.com/smarco/gem3-mapper.git -b ${GEM_MAPPER_VERSION} ./ && \
    chmod +x configure && \
    ./configure && \
    make && \
    mv ${SRC_BASE}/bin ${INSTALL_BASE} && \
    ln -s ${INSTALL_BASE}/bin/* /usr/local/bin/

RUN which gem-mapper 

