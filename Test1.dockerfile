FROM rocker/r-ver:4.3.1
RUN apt-get update && apt-get install -y libglpk40 libgsl-dev zlib1g-dev libcairo2-dev libx11-dev libxt-dev liblzma-dev libbz2-dev libcurl4-openssl-dev wget && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages(c('Seurat', 'scales', 'cowplot', 'ggplot2', 'BiocManager', 'dplyr', 'viridis', 'grid', 'Libra', 'rlang', 'patchwork', 'reshape2', 'magrittr', 'tibble', 'purrr', 'forcats', 'remotes'))"
RUN R -e "BiocManager::install(version = '3.17')"
RUN R -e "BiocManager::install(c('biomaRt', 'edgeR', 'SingleCellExperiment', 'scDblFinder'))"
WORKDIR /opt
RUN wget -O cellranger-arc-2.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-arc/cellranger-arc-2.0.2.tar.gz?Expires=1695451161&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hcmMvY2VsbHJhbmdlci1hcmMtMi4wLjIudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjk1NDUxMTYxfX19XX0_&Signature=hWypVLHLgNXR4OvZ5GZGAAc6xpkkBzXYY7Irb1WJi8gC0Fwwt9k6-h-WUfyrHbYS6Wcxmz~tKmfel9PBAeEHwUU3sh3myYQq4prtL4pKvRTfcySFT9IiHMZbE1hZW-9XCrTWmbYNPTreT0sGooRaActBwJpNpVyiOTkxx8xHQnEGHoCbBTrptGjehJIsOoDsgMZESFydecO2Zdr1vHFKpdPR-DG2wBuDHt0WKO0lstZHzlK1J0~A7ze-9q-0rTWQP0YR0fEg8i3tmAYxf-ueCwI2NnABvsVZtqUrlw5cTB9H27IylQIHJoy8zfEBrlS0g7LuZ7KccIRS9GIeTR3tag__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
RUN tar -xzvf cellranger-arc-2.0.2.tar.gz && rm cellranger-arc-2.0.2.tar.gz
ENV PATH="/opt/cellranger-arc-2.0.2:$PATH"
RUN apt-get update && apt-get install -y git libfftw3-dev && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o /bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member
RUN apt-get update && apt-get install -y python3.11 pip && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN pip install cellbender
CMD ["R"]
