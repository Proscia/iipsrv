#
# Copyright (c) 2009-2018. Authors: see NOTICE file.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

FROM ubuntu:18.04

SHELL ["/bin/bash", "-c"]

RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list

RUN cat /etc/apt/sources.list

RUN apt-get -y update && \
    apt-get install -y \
      cron \
      git \
      libjpeg8-dev \
      libmemcached-dev \
      libopenjp2-7-dev \
      libssl-dev \
      libtiff5-dev \
      libpng-dev \
      psmisc \
      autoconf \
      libtool \
      rsyslog \
      spawn-fcgi \
      cmake \
      wget && \
    apt-get install -y \
      software-properties-common && \
    apt-get -y build-dep openslide

ARG NGINX_VERSION=1.15.1
RUN cd /tmp && \
    wget https://nginx.org/download/nginx-$NGINX_VERSION.tar.gz && \
    tar -zxvf nginx-$NGINX_VERSION.tar.gz && \
    rm nginx-$NGINX_VERSION.tar.gz && \
    cd /tmp/nginx-$NGINX_VERSION && \
    ./configure && \
    make && \
    make install && \
    rm -rf /tmp/*

# openslide
ARG OPENSLIDE_VERSION=3.4.1
RUN cd /tmp && \
    wget https://github.com/openslide/openslide/releases/download/v${OPENSLIDE_VERSION}/openslide-${OPENSLIDE_VERSION}.tar.gz && \
    tar -zxvf ./openslide-${OPENSLIDE_VERSION}.tar.gz && \
    cd ./openslide-${OPENSLIDE_VERSION} && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install && \
    ldconfig

RUN cd /opt && \
    git clone https://5af320861a6c1ff842a80918fc4e846b8ffa0fae:x-oauth-basic@github.com/Proscia/libCZI.git && \
    cd libCZI && \
    cmake -G "Unix Makefiles" && \
    make && \
    shopt -s globstar && \
    ls -l **/*.{a,so} **/CZIcmd && \
    [ -f Src/libCZI/liblibCZI.so ] \
        && [ -f Src/libCZI/liblibCZIStatic.a ] \
        && [ -f Src/JxrDecode/libJxrDecodeStatic.a ] \
        && [ -x Src/CZICmd/CZIcmd ]

#IIP Image Server
RUN cd /opt && \
    git clone https://5af320861a6c1ff842a80918fc4e846b8ffa0fae:x-oauth-basic@github.com/Proscia/iipsrv.git
RUN cd /opt/iipsrv && \
    git checkout optimize-plus-docker
RUN cd /opt/iipsrv && \
    sh autogen.sh && \
    ./configure && \
    make

COPY deploy.sh /tmp/deploy.sh
RUN chmod +x /tmp/deploy.sh

COPY start-iip.sh /tmp/start-iip.sh
RUN chmod +x /tmp/start-iip.sh

COPY check-status.sh /tmp/check-status.sh
RUN chmod +x /tmp/check-status.sh

RUN ldconfig -v

ENTRYPOINT ["/tmp/deploy.sh"]
