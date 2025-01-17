# Instructions for building a container for testing
# or experimenting with the Tango end-to-end simulator. The main 


# This container sets up the environment by installing many scientific
# libraries that SPEXone depends on. That said, it doesn't have to
# limited to SPEXone. The list can easily be extended by editing
# spack.yaml, which lists packages installed by the Spack package
# manager.

FROM ubuntu:24.04 AS builder-ubuntu

ENV DEBIAN_FRONTEND=noninteractive
RUN export DEBCONF_NOWARNINGS=yes \
    && apt-get -yq update \
    && apt-get -yq upgrade \
    && apt-get install -y apt-utils locales \
    && sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen \
    && locale-gen \
    && apt-get -yq install --no-install-recommends \
    build-essential \
    ca-certificates \
    clang-format-18 \
    curl \
    file \
    g++ \
    gcc \
    gfortran \
    git \
    gnupg2 \
    iproute2 \
    make \
    python3 \
    python3-pip \
    python3-setuptools \
    unzip \
    wget \
    liblapack-dev \
    libeigen3-dev \
    libhdf5-dev
ENV LANGUAGE=en_US.UTF-8 LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

FROM builder-ubuntu AS builder-science

# Set up Spack
ENV spack_version=0.22.1
WORKDIR /opt
RUN DL_URL=https://github.com/spack/spack/releases/download/ \
    && wget $DL_URL/v$spack_version/spack-$spack_version.tar.gz \
    && tar xf spack-$spack_version.tar.gz \
    && . spack-$spack_version/share/spack/setup-env.sh \
    && spack spec curl

# Build Spack environment
COPY spack.yaml ./spack-environment/spack.yaml
RUN . ./spack-$spack_version/share/spack/setup-env.sh \
    && spack env activate spack-environment \
    && spack install --fail-fast \
    && spack gc -y

# Build TEDS Python requirements
RUN git clone https://bitbucket.org/sron_earth/teds.git \
    && cd teds \
    && git checkout 32f3847bfc \
    && /opt/view/bin/python -m venv /venv \
    && . /venv/bin/activate \
    && pip install --upgrade pip \
    && pip install -e .

# Expose the environment
RUN . ./spack-$spack_version/share/spack/setup-env.sh \
    && spack env activate --sh -d spack-environment >> /etc/profile.d/z10_spack_environment.sh
SHELL ["/bin/bash", "--rcfile", "/etc/profile", "-l", "-c"]

FROM builder-ubuntu

COPY --from=builder-science /opt/software /opt/software
COPY --from=builder-science /opt/._view /opt/._view
COPY --from=builder-science /opt/view /opt/view
COPY --from=builder-science /etc/profile.d /etc/profile.d
COPY --from=builder-science /venv /venv

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l", "-c"]
CMD ["/bin/bash"]
