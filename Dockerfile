FROM nfcore/base
LABEL authors="QBiC" \
      description="Docker image containing all requirements for nf-core/rnadeseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-rnadeseq-1.0dev/bin:$PATH
