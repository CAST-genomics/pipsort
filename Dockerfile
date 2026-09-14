FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
	build-essential \
	cmake \
	libgsl-dev \
	libopenblas-dev \
	git-all

RUN git clone https://github.com/cast-genomics/pipsort
WORKDIR pipsort/
RUN git checkout cmake
RUN mkdir build
WORKDIR build/
RUN cmake ..
RUN make
RUN cmake --install . --prefix /usr
