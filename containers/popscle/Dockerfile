# popscle built from statgen/popscle master @ da70fc78 — the EXACT commit the
# retired Demuxafy.sif shipped (and that the brca_mega_atlas results used).
#
# Why from source (not a biocontainer): the only public popscle biocontainer is
# the older `0.1beta` release, whose `dsc-pileup` emits the per-cell UMI table in
# a different format (`chr:start-end`) than master @ da70fc78 (`chr:pos:count:+`).
# Matching the commit keeps pileup.umi.gz byte-identical to the atlas golden.
# bioconda has no package for this commit, so this is the one tool that needs a
# Dockerfile (cf. imputation's containers/minimac4). The other engines use public
# biocontainers / the souporcell image — see config.yaml `containers:`.
#
# Mirrors Demuxafy's bundled /opt/popscle/Dockerfile (statgen/popscle + a
# from-source htslib), pinned for reproducibility. ubuntu:20.04 (gcc 9) is
# deliberate: popscle @ da70fc78 omits an `#include <limits>` that gcc >= 11
# (ubuntu 22.04) rejects ("'numeric_limits' is not a member of 'std'"); gcc 9
# includes it transitively and builds the commit as Demuxafy's toolchain did.
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ARG HTSLIB_TAG=1.16
ARG POPSCLE_COMMIT=da70fc78da385ef049e0e890342acfd62842cae0

RUN apt-get update && apt-get install -y --no-install-recommends \
        autoconf automake build-essential cmake git \
        libbz2-dev libcurl4-openssl-dev liblzma-dev libssl-dev zlib1g-dev \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# htslib (popscle links against it; built from a pinned release tag)
RUN git clone --depth 1 --branch "${HTSLIB_TAG}" --recurse-submodules \
        https://github.com/samtools/htslib.git /tmp/htslib \
 && cd /tmp/htslib \
 && autoreconf -i && ./configure --prefix=/usr/local && make -j"$(nproc)" && make install \
 && ldconfig && rm -rf /tmp/htslib

# popscle pinned to the Demuxafy commit
RUN git clone https://github.com/statgen/popscle.git /tmp/popscle \
 && cd /tmp/popscle \
 && git checkout "${POPSCLE_COMMIT}" \
 && mkdir build && cd build && cmake .. && make -j"$(nproc)" \
 && cp /tmp/popscle/bin/popscle /usr/local/bin/popscle \
 && chmod +x /usr/local/bin/popscle \
 && rm -rf /tmp/popscle

CMD ["popscle"]
