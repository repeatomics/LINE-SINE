FROM continuumio/miniconda3:latest

# system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install bioinformatics tools
RUN conda install -c bioconda -c conda-forge -y \
    bwa=0.7.17 \
    samtools=1.17 \
    bcftools=1.18 \
    fastp=0.23.4 \
    python=3.10 \
    && conda clean -a -y

# python dependencies
RUN pip install --no-cache-dir \
    pandas==2.0.3 \
    numpy==1.24.3 \
    matplotlib==3.7.2 \
    seaborn==0.12.2 \
    scikit-learn==1.3.0 \
    pysam==0.21.0 \
    jupyter

# working directory
WORKDIR /app

# copy project
COPY . /app/

# entrypoint
ENTRYPOINT ["python", "pipeline.py"]
