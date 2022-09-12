FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get install -y cmake ninja-build xxd clang python3.8

RUN mkdir host

COPY cg cg
COPY CMakeLists.txt CMakeLists.txt
COPY build.py build.py

RUN python3 build.py build

WORKDIR "/host"
ENTRYPOINT ["/build/cg/cg"]
