FROM ubuntu:bionic

#RUN /bin/sh -c [ -z "$(apt-get indextargets)" ]

RUN set -xe 		&& echo '#!/bin/sh' > /usr/sbin/policy-rc.d 	&& echo 'exit 101' >> /usr/sbin/policy-rc.d 	&& chmod +x /usr/sbin/policy-rc.d 		&& dpkg-divert --local --rename --add /sbin/initctl 	&& cp -a /usr/sbin/policy-rc.d /sbin/initctl 	&& sed -i 's/^exit.*/exit 0/' /sbin/initctl 		&& echo 'force-unsafe-io' > /etc/dpkg/dpkg.cfg.d/docker-apt-speedup 		&& echo 'DPkg::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' > /etc/apt/apt.conf.d/docker-clean 	&& echo 'APT::Update::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' >> /etc/apt/apt.conf.d/docker-clean 	&& echo 'Dir::Cache::pkgcache ""; Dir::Cache::srcpkgcache "";' >> /etc/apt/apt.conf.d/docker-clean 		&& echo 'Acquire::Languages "none";' > /etc/apt/apt.conf.d/docker-no-languages 		&& echo 'Acquire::GzipIndexes "true"; Acquire::CompressionTypes::Order:: "gz";' > /etc/apt/apt.conf.d/docker-gzip-indexes 		&& echo 'Apt::AutoRemove::SuggestsImportant "false";' > /etc/apt/apt.conf.d/docker-autoremove-suggests

RUN mkdir -p /run/systemd && echo 'docker' > /run/systemd/container

RUN apt-get update \
 && apt-get install -y \
 software-properties-common \
 build-essential \
 && apt-get update

ENV DEBIAN_FRONTEND="noninteractive" TZ="Asia/Tokyo"

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get update \
  && apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  && mkdir /root/.conda \
  && bash Miniconda3-latest-Linux-x86_64.sh -b \
  && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

RUN apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y --fix-missing \
 git \
 nano \
 r-base r-recommended \
 samtools \
 htop \
 pigz \
 wget \
 curl \
 libcurl4-openssl-dev \
 libcairo2-dev cairo-dock cairo-dock-plug-ins \
 libxt-dev \
 libharfbuzz-dev libfribidi-dev \
 libhdf5-dev \
 libgl1-mesa-glx libegl1-mesa \
 libxrandr2 libxrandr2 \
 libxss1 libxcursor1 libxcomposite1 \
 libasound2 libxi6 libxtst6 \
 libssl-dev \
 python-pip

RUN conda install numpy scipy cython numba matplotlib scikit-learn h5py click

RUN pip install pysam velocyto

RUN apt update -qq \
  && apt install --no-install-recommends software-properties-common dirmngr \
  && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
  && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/" \
  && apt install -y --no-install-recommends r-base

RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' \
  && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
  && apt update \
  && gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
  && gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | apt-key add - \
  && apt-get install -y r-base r-base-core r-recommended r-base-dev

RUN R -e 'install.packages(c("yaml", "BiocManager", "Cairo", "ggplot2", "cowplot", "data.table", "devtools", "dplyr", "ggraster", "ggrastr", "hdf5r", "inflection", "ragg", "textshaping"), dependencies = TRUE)'

RUN R -e 'BiocManager::install(c("GenomicFeatures", "plyranges", "Rhtslib", "Rsamtools"), ask = FALSE)'

RUN R -e 'install.packages("devtools"); devtools::install_github(repo = "hhoeflin/hdf5r"); devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")'

RUN mkdir -p home && cd home && git clone https://github.com/TomKellyGenetics/zUMIs.git

RUN cd /home  && wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz && tar -xzf 2.7.9a.tar.gz && cd STAR-2.7.9a/source && make

CMD ["/bin/bash"]

#bash
