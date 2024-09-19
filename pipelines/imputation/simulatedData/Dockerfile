FROM python:3.12

# update and install gsl
RUN apt-get -y update && apt-get -y install libgsl-dev  && apt-get install -y tabix

# install requirements for stdpopsim and stdpopsim
RUN pip install msprime --no-binary msprime && pip install scipy && pip install stdpopsim && pip install tskit
