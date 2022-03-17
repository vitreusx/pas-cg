FROM ubuntu:20.04

RUN apt-get update && apt-get install -y cmake ninja-build xxd g++

ARG BUILD_TYPE
RUN test -n "${BUILD_TYPE}" \
    || { echo "BUILD_TYPE must be set!"; exit 1; }

COPY . .

RUN mkdir build \
    && cd build \
    && cmake -DCMAKE_C_COMPILER=/usr/bin/gcc \
             -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
             -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
             -G Ninja .. \
    && ninja cg

