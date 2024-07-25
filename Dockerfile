### Copyright 2003-2024. GenePattern Team @ Mesirov Lab - University of California, San Diego. All rights reserved.
#
# Currently, module uses python:3.11 image.
FROM python:3.11

# Based off of ExampleModule Dockerfile found at https://github.com/genepattern/ExampleModule 
LABEL maintainer="Julia Kononova jkononova@ucsd.edu, Omar Halawa ohalawa@ucsd.edu"

# Setting up proper environment, see ExampleModule Dockerfile for more info
# -----------------------------------
# -----------------------------------

# Ensuring up-to-date pip and importing necessary modules 
RUN pip install --upgrade pip && \
    pip install numpy==1.23.5 \
     pandas==2.2.0 \
     scipy==1.11.1 \
     jupyterlab==3.6.3 \
     matplotlib==3.7.2 \
     scikit-learn==1.3.0 \
     scanpy==1.9.8 \
     tqdm==4.65.0 \
     seaborn==0.13.2 \
     squidpy==1.4.1 \
     tangram-sc==1.0.4
     
# Build using "docker build -t <TAG> ."
# Run using "docker run -it --rm <IMAGE ID> bash"

RUN useradd -ms /bin/bash gpuser
USER gpuser
WORKDIR /home/gpuser

USER root
RUN mkdir /Tangram \
    && chown gpuser /Tangram

USER gpuser
COPY src/*.py /Tangram/