# hash:sha256:481e8cc42c19dda7ed261a03cf28df8633483d35ebe4da8d6ba1c9f382cd1609
FROM registry.codeocean.com/codeocean/gcc:5.4.0-ubuntu16.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gfortran=4:5.3.1-1ubuntu1 \
    && rm -rf /var/lib/apt/lists/*
