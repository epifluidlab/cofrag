# A minimal Docker container for cofrag project

FROM continuumio/miniconda3:4.8.2

RUN conda config --add channels bioconda && conda config --add channels conda-forge
RUN conda install -y r-base=3.6.3
RUN conda install -y r-essentials=3.6
RUN conda install -y r-optparse=1.6.6
RUN conda install -y r-here=0.1
RUN conda install -y r-doparallel=1.0.15
RUN conda install -y r-logging=0.10
RUN conda install -y -c bioconda bioconductor-genomicranges=1.38.0