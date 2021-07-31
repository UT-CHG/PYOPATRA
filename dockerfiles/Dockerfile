FROM ubuntu:20.04
RUN apt-get update
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install libmpich12 \
    libmpich-dev \
    mpich \
    libhdf5-dev \
    cmake \
    libeigen3-dev \
    pybind11-dev \
    python3 \
    libpython3-dev \
    python3-pip \
    python-is-python3 \
    wget \
    && wget https://github.com/catchorg/Catch2/archive/v2.13.4.tar.gz \
    && tar -xzf v2.13.4.tar.gz \
    && cd Catch2-2.13.4 \
    && cmake -Bbuild -H. -DBUILD_TESTING=OFF \
    && cmake --build build/ --target install \
    && cd .. \
    && pip3 install numpy scipy matplotlib coverage pytest netcdf4 mpi4py \
    && mkdir data
COPY . /PYOPATRA
RUN pip3 install PYOPATRA/
