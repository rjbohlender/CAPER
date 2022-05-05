FROM debian:bullseye AS builder
# FROM alpine AS builder

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && apt-get install -y --no-install-recommends apt-utils
# RUN apk update

RUN apt-get install -y --no-install-recommends \
     ca-certificates \
     file \
     sudo \
     unzip \
     wget
RUN apt-get install -y --no-install-recommends build-essential cmake libboost-all-dev libarmadillo-dev
# RUN apk add make g++ cmake

# RUN apk add boost-dev armadillo-dev blas-dev lapack-dev

COPY . /src/

WORKDIR /src/build

RUN cmake -DCMAKE_BUILD_TYPE=Release ../.

RUN cmake --build . --target carva

FROM debian:bullseye

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /

COPY --from=builder /src/build/carva /bin
RUN mkdir -p /filter
COPY --from=builder /src/filter/filter_whitelist.csv /filter/filter_whitelist.csv

RUN apt-get update

RUN apt-get install -y --no-install-recommends libboost-program-options1.74.0 libboost-iostreams1.74.0 libarmadillo10
