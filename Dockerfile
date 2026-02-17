FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc-14 \
        libc6-dev \
        make \
        libfftw3-dev \
        libsndfile1-dev \
        portaudio19-dev \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY . .
CMD ["sh", "-c", "make clean && make && make test"]
