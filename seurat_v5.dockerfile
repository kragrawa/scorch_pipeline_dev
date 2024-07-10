FROM rocker/r-ver:4.3.1
RUN apt-get update && apt-get install -y libglpk40 libgsl-dev zlib1g-dev libcairo2-dev libx11-dev libxt-dev liblzma-dev libbz2-dev libcurl4-openssl-dev wget && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages(c('scales', 'cowplot', 'ggplot2', 'BiocManager', 'dplyr', 'viridis', 'grid', 'Libra', 'rlang', 'patchwork', 'reshape2', 'magrittr', 'tibble', 'purrr', 'forcats', 'remotes'))"
RUN R -e "remotes::install_github('satijalab/seurat', 'seurat5', quiet = TRUE)"
RUN R -e "remotes::install_github('satijalab/seurat-data', 'seurat5', quiet = TRUE)"
RUN R -e "remotes::install_github('satijalab/azimuth', quiet = TRUE)"
RUN R -e "remotes::install_github('satijalab/seurat-wrappers', quiet = TRUE)"
RUN R -e "remotes::install_github('stuart-lab/signac', quiet = TRUE)"
RUN R -e "BiocManager::install(version = '3.18')"
RUN R -e "BiocManager::install(c('biomaRt', 'edgeR', 'SingleCellExperiment', 'scDblFinder'))"
WORKDIR /opt
#RUN wget -O cellranger-arc-2.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-arc/cellranger-arc-2.0.2.tar.gz?Expires=1708576779&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hcmMvY2VsbHJhbmdlci1hcmMtMi4wLjIudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNzA4NTc2Nzc5fX19XX0_&Signature=PNQAarAKzP4hxFGr~giT~duX7kebKvKg4X28RUBK~d4ZmQURkMDRTyyMSnWZ1c~xUak2WQKnvUptM1h8xZQ3G9YF7xFPlQ75tisnYnTwkifsPlNsMVka~eOJr1wJmok3hPkaQPxqcayjmW8rrZ9umZnwvYCsIPtsE98DeDoAjlW1QKwhaRZDyGT5FyucLwgG7UGPHsXtbWAjtkRh8mAK3zd10Ief6vbHmAj6AWXMHJ~Atp8gAnS4T~BpYMj0IfbNLxXnaynQ9Auz2DmcedtB0~YHSuHucKoKeOW1o3Q3ZcV-TGrPnrHC~21726AUV5pPHuYT2JoM7A2AjK1G7tPIrA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
#RUN tar -xzvf cellranger-arc-2.0.2.tar.gz && rm cellranger-arc-2.0.2.tar.gz
#ENV PATH="/opt/cellranger-arc-2.0.2:$PATH"
RUN apt-get update && apt-get install -y git libfftw3-dev && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o /bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member
RUN apt-get update && apt-get install -y python3.11 pip && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN pip install cellbender
CMD ["R"]