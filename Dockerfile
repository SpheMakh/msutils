FROM kernsuite/base:1
RUN docker-apt-install python-casacore \ 
    python-numpy python-owlcat \
    meqtrees-timba

ADD . /src
RUN pip install /src

