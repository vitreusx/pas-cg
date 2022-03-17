FROM ubuntu:20.04

RUN apt-get update && apt-get install -y cmake ninja-build xxd g++

RUN mkdir host

COPY . .

RUN mkdir build \
    && cd build \
    && cmake -DCMAKE_C_COMPILER=/usr/bin/gcc \
             -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
             -DCMAKE_BUILD_TYPE=Release \
             -G Ninja .. \
    && ninja cg

WORKDIR "/host"
ENTRYPOINT ["/build/cg/cg"]