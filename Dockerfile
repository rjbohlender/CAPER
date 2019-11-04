FROM debian:bullseye AS builder

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get install -y --no-install-recommends \
      ca-certificates \
      file \
      sudo \
      unzip \
      wget

RUN apt-get install -y --no-install-recommends \
      build-essential \
      cmake \
      lcov

RUN apt-get install -y --no-install-recommends libboost-program-options-dev libboost-iostreams-dev libarmadillo-dev cmake

COPY . /src/

WORKDIR /src/build

RUN cmake -DCMAKE_BUILD_TYPE=Release ../.

RUN make carva
RUN make vcf2matrix

FROM debian:bullseye

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /

COPY --from=builder /src/build/carva /bin

RUN apt-get -y update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get install -y --no-install-recommends libboost-program-options-dev libboost-iostreams-dev libarmadillo-dev
