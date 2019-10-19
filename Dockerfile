FROM debian:bullseye

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get install -y --no-install-recommends \
      bash-completion \
      bzip2 \
      ca-certificates \
      curl \
      file \
      git \
      patch \
      procps \
      sudo \
      unzip \
      vim \
      wget

RUN apt-get install -y --no-install-recommends \
      build-essential \
      clang++-6 \
      cmake \
      gdb \
      lcov \
      valgrind

# asio
RUN apt-get install -y --no-install-recommends \
      libssl-dev

# beast
RUN apt-get install -y --no-install-recommends \
      libssl-dev \
      zlib1g-dev

# gil
RUN apt-get install -y --no-install-recommends \
      libjpeg-dev \
      libpng-dev \
      libtiff-dev \
      zlib1g-dev

# integer
RUN apt-get install -y --no-install-recommends \
      libgmp-dev

# iostreams
RUN apt-get install -y --no-install-recommends \
      libbz2-dev \
      liblzma-dev \
      zlib1g-dev

# locale
RUN apt-get install -y --no-install-recommends \
      libicu-dev

# mpi
RUN apt-get install -y --no-install-recommends \
      libopenmpi-dev

# python
RUN apt-get install -y --no-install-recommends \
      python3-dev

RUN apt-get install -y --no-install-recommends libboost-all-dev libarmadillo-dev cmake

COPY . /src/

WORKDIR /src/build

RUN cmake ../.

RUN make
