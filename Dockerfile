# FROM debian:bullseye AS builder
FROM alpine AS builder

# ENV DEBIAN_FRONTEND noninteractive

# RUN apt-get -y update && apt-get install -y --no-install-recommends apt-utils
RUN apk update

# RUN apk add \
#      ca-certificates \
#      file \
#      sudo \
#      unzip \
#      wget

RUN apk add make g++ cmake

RUN apk add boost-dev armadillo-dev blas-dev lapack-dev

COPY . /src/

WORKDIR /src/build

RUN cmake -DCMAKE_BUILD_TYPE=Release ../.

RUN make carva

FROM alpine

# ENV DEBIAN_FRONTEND noninteractive

WORKDIR /

COPY --from=builder /src/build/carva /bin

RUN apk update

RUN apk add boost blas lapack armadillo
