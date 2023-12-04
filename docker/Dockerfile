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

FROM ubuntu:22.04

SHELL ["/bin/bash", "-c"]

ARG QUIET="perl -ne '\$|=1; print (++\$n%80==1 ? \"\\n\$_\" : \".\");'"

RUN echo "QUIET: $QUIET" && \
    sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list && \
    cat /etc/apt/sources.list

RUN apt-get update -y \
      | eval "$QUIET" && \
    apt-get upgrade -y | eval "$QUIET" && \
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
      wget \
      | eval "$QUIET" && \
    apt-get install -y \
      software-properties-common \
      | eval "$QUIET" && \
    apt-get build-dep -y openslide \
      | eval "$QUIET"


# openslide
ARG OPENSLIDE_VERSION=3.4.1
RUN cd /tmp && \
    git clone https://github.com/Proscia/openslide.git && \
    cd ./openslide && \
    autoreconf -i && \
    ./configure && \
    make -j$(nproc) && \
    make install && \
    ldconfig


# Zeiss CZI
ARG LIBCZI_REPOSITORY_URL=https://5af320861a6c1ff842a80918fc4e846b8ffa0fae:x-oauth-basic@github.com/Proscia/libCZI.git
ARG LIBCZI_BRANCH=master

RUN cd /opt && git clone ${LIBCZI_REPOSITORY_URL} && \
    cd /opt/libCZI && git checkout ${LIBCZI_BRANCH} && git remote -v && git branch

RUN cd /opt/libCZI && \
    cmake -G "Unix Makefiles" && \
    make && \
    shopt -s globstar && ls -l **/*.{a,so} **/CZIcmd && \
    [ -f Src/libCZI/liblibCZI.so ] \
        && [ -f Src/libCZI/liblibCZIStatic.a ] \
        && [ -f Src/JxrDecode/libJxrDecodeStatic.a ] \
        && [ -x Src/CZICmd/CZIcmd ]

COPY . /opt/iipsrv
RUN cd /opt/iipsrv && ls

# RUN cd /opt && git clone ${IIPSRV_REPOSITORY_URL} && \
#     cd /opt/iipsrv && git checkout ${IIPSRV_BRANCH} && git remote -v && git branch

RUN cd /opt/iipsrv && \
    sh autogen.sh && \
    ./configure && \
    make && \
    shopt -s globstar && ls -l **/*.fcgi


COPY docker/deploy.sh /tmp/deploy.sh
RUN chmod +x /tmp/deploy.sh

COPY docker/start-iip.sh /tmp/start-iip.sh
RUN chmod +x /tmp/start-iip.sh

COPY docker/check-status.sh /tmp/check-status.sh
RUN chmod +x /tmp/check-status.sh

RUN ldconfig -v

RUN useradd worker --uid 1000
RUN chown worker:worker /tmp
USER worker:worker

ENTRYPOINT ["/tmp/deploy.sh"]
