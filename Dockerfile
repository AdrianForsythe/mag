FROM nfcore/base
=======
MAINTAINER Hadrien Gourlé <hadrien.gourle@slu.se>
LABEL description="Docker image containing all requirements for nf-core/mag pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-mag-1.0dev/bin:$PATH
