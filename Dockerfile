FROM condaforge/mambaforge
LABEL org.opencontainers.image.source=https://github.com/qbic-pipelines/rnadeseq
LABEL org.opencontainers.image.description="Docker image containing all requirements for qbic-pipelines/rnadeseq pipeline"
LABEL org.opencontainers.image.authors="Gisela Gabernet, Alexander Peltzer, Oskar Wacker"
LABEL org.opencontainers.image.licenses=MIT
COPY environment.yml /
#RUN conda install -c conda-forge mamba
RUN mamba env create --file /environment.yml -p /opt/conda/envs/qbic-pipelines-rnadeseq-2.5 && \
    mamba clean --all --yes
RUN apt-get update -qq && \
    apt-get install -y zip procps ghostscript
# Add conda installation dir to PATH
ENV PATH /opt/conda/envs/qbic-pipelines-rnadeseq-2.5/bin:$PATH
# Dump the details of the installed packates to a file for posterity
RUN mamba env export --name qbic-pipelines-rnadeseq-2.5 > qbic-pipelines-rnadeseq-2.5.yml
# Instruct R processes to use these empty files instead of clashing with a local config
RUN touch .Rprofile
RUN touch .Renviron
