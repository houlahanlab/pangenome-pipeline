FROM ubuntu:16.04 AS builder

RUN apt-get -yq update \
    && apt-get install -yq --no-install-recommends \
    locales \
    ca-certificates \
    gcc \
    g++ \
    bzip2 \
    zlib1g-dev \
    make \
    cmake \
    python \
    python3 \
    wget \
    unzip \
    tabix \
    libboost-dev libboost-date-time-dev libboost-filesystem-dev libboost-program-options-dev libboost-serialization-dev \
    libboost-system-dev libboost-timer-dev libboost-chrono-dev libboost-test-dev \
    autoconf automake perl libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT=/opt/tools
ENV PERL5LIB=$OPT/lib/perl5
ENV LD_LIBRARY_PATH=$OPT/lib
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

# download and install htslib
ENV SAMTOOLS_VERSION="1.23.1"

WORKDIR /tmp
RUN wget --progress=dot:giga https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -jxf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd /tmp/htslib-${SAMTOOLS_VERSION}/ \
    && ./configure --without-libdeflate --prefix=$OPT/htslib \
    && make \
    && make install

# download and install samtools
WORKDIR /tmp
RUN wget --progress=dot:giga https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -jxf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd /tmp/samtools-${SAMTOOLS_VERSION} \
    && ./configure --without-libdeflate --prefix=$OPT/samtools \
    && make \
    && make install

# download and install bcftools
WORKDIR /tmp
RUN wget --progress=dot:giga https://github.com/samtools/bcftools/releases/download/${SAMTOOLS_VERSION}/bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -jxf bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd /tmp/bcftools-${SAMTOOLS_VERSION} \
    && ./configure --without-libdeflate --prefix=$OPT/bcftools \
    && make \
    && make install

# clean up
WORKDIR /tmp
RUN rm -rf htslib-${SAMTOOLS_VERSION}* \
    && rm -rf samtools-${SAMTOOLS_VERSION}* \
    && rm -rf bcftools-${BCFTOOLS_VERSION}*

# download and install strelka2
ENV STRELKA_VERSION="2.9.10"

WORKDIR /tmp
RUN wget --progress=dot:giga https://github.com/Illumina/strelka/releases/download/v${STRELKA_VERSION}/strelka-${STRELKA_VERSION}.centos6_x86_64.tar.bz2 \
    && tar -jxf strelka-${STRELKA_VERSION}.centos6_x86_64.tar.bz2 \
    && mkdir $OPT/strelka \
    && cp -r strelka-${STRELKA_VERSION}.centos6_x86_64/* $OPT/strelka \
    && rm -rf strelka-${STRELKA_VERSION}*

# download and install gatk
ENV GATK_VERSION="4.3.0.0"

WORKDIR /tmp
RUN wget --progress=dot:giga https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
    && unzip gatk-${GATK_VERSION}.zip \
    && mkdir -p $OPT/gatk \
    && cp gatk-${GATK_VERSION}/gatk-package-4.3.0.0-local.jar $OPT/gatk/ \
    && cp gatk-${GATK_VERSION}/gatk $OPT/gatk/ \
    && rm gatk-${GATK_VERSION}.zip

# download and install vg
ENV VG_VERSION="1.61.0"

WORKDIR $OPT/
RUN wget --progress=dot:giga https://github.com/vgteam/vg/releases/download/v${VG_VERSION}/vg && \
    chmod +x vg

FROM ubuntu:16.04

LABEL maintainer="prokos2@mcmaster.ca" \
      version="1.0.0" \
      description="Pangenome:Somatic pipeline"

RUN apt-get -yq update \
    && apt-get install -yq --no-install-recommends \
    locales \
    bzip2 \
    zlib1g \
    python \
    tabix \
    libcurl4-gnutls-dev \
    openjdk-8-jre-headless \
    unattended-upgrades && \
    unattended-upgrade -d -v && \
    apt-get remove -yq unattended-upgrades && \
    apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT=/opt/tools
ENV PERL5LIB=$OPT/lib/perl5
ENV LD_LIBRARY_PATH=$OPT/lib
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

ENV PATH=$OPT/samtools/bin:$PATH
ENV PATH=$OPT/bcftools/bin:$PATH
ENV PATH=$OPT/strelka/bin:$PATH
ENV PATH=$OPT/gatk:$PATH
ENV PATH=$OPT:$PATH

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]

# test to verify tool binaries work
RUN configureStrelkaSomaticWorkflow.py --version && \
    samtools --version && \
    bcftools --version && \
    gatk --version

#    vg version

