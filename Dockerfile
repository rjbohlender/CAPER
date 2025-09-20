FROM debian:bullseye AS builder
# FROM alpine AS builder

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        apt-utils \
        ca-certificates \
        file \
        sudo \
        unzip \
        wget \
        build-essential \
        cmake \
        libboost-all-dev \
        libarmadillo-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# RUN apk update
# RUN apk add make g++ cmake

# RUN apk add boost-dev armadillo-dev blas-dev lapack-dev

COPY . /src/

WORKDIR /src/build

RUN cmake -DCMAKE_BUILD_TYPE=Release ../.

RUN cmake --build . --target caper

FROM debian:bullseye

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /

COPY --from=builder /src/build/caper /bin
RUN mkdir -p /filter
COPY --from=builder /src/filter/filter_whitelist.csv /filter/filter_whitelist.csv

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libboost-program-options1.74.0 \
        libboost-iostreams1.74.0 \
        libarmadillo10 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
