FROM kernsuite/base:5
RUN docker-apt-install python3-casacore python3-pip
ADD . /msutils-src
RUN pip3 install /msutils-src
