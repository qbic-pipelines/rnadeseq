FROM nfcore/base
LABEL authors="QBiC" \
      description="Docker image containing all requirements for qbicsoftware/rnadeseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN apt-get install -y zip
ENV PATH /opt/conda/envs/nf-core-rnadeseq-1.0dev/bin:$PATH
