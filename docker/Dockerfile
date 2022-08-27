FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get install -y cmake ninja-build xxd g++

RUN mkdir host

COPY cg cg
COPY CMakeLists.txt CMakeLists.txt

RUN mkdir build \
    && cd build \
    && cmake -DCMAKE_C_COMPILER=/usr/bin/gcc \
             -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
             -DCMAKE_BUILD_TYPE=Release \
             -G Ninja .. \
    && ninja cg

WORKDIR "/host"
ENTRYPOINT ["/build/cg/cg"]
