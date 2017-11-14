FROM kernsuite/base:3
RUN docker-apt-install python-casacore \ 
    meqtrees
RUN docker-apt-install python-pip
RUN pip install -U pip
ADD . /msutils-src
RUN pip install /msutils-src
