FROM condaforge/mambaforge
LABEL authors="Gisela Gabernet, Alexander Peltzer" \
    description="Docker image containing all requirements for qbic-pipelines/rnadeseq pipeline"
COPY environment.yml /
#RUN conda install -c conda-forge mamba
RUN mamba env create --file /environment.yml -p /opt/conda/envs/qbic-pipelines-rnadeseq-2.1 && \
    mamba clean --all --yes
RUN apt-get update -qq && \
    apt-get install -y zip procps ghostscript
# Add conda installation dir to PATH
ENV PATH /opt/conda/envs/qbic-pipelines-rnadeseq-2.1/bin:$PATH
# Dump the details of the installed packates to a file for posterity
RUN mamba env export --name qbic-pipelines-rnadeseq-2.1 > qbic-pipelines-rnadeseq-2.1.yml
# Instruct R processes to use these empty files instead of clashing with a local config
RUN touch .Rprofile
RUN touch .Renviron

