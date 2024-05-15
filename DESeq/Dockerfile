# Use the official R base image
FROM rocker/r-ver:4.0.3

RUN apt-get update && apt-get install -y \
    build-essential \
    pandoc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('dplyr','knitr','gridextra', 'tidyverse', 'data.table', 'foreach', 'rmarkdown'), repos = 'https://cloud.r-project.org/')"

# add bioconductor packages not in cran
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

ADD rpackages.R /tmp/
RUN R -f /tmp/rpackages.R && rm -r /tmp/rpackages.R


# Set the working directory
WORKDIR /data

# Copy the local code to the working directory in the Docker container
COPY /tx2gene_hg38.csv /data
# RUN chmod +x /data/pre_processing.py

# TODO: add my scripts ?