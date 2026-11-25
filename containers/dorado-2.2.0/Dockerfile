# syntax=docker/dockerfile:1.7@sha256:a57df69d0ea827fb7266491f2813635de6f17269be881f696fbfdf2d83dda33e

FROM nvidia/cuda:12.8.1-cudnn-runtime-ubuntu24.04@sha256:9175fa92f96de35a8cfb9493f0dfcf9435c7a597e9d95ad41d2cae382a95e3f9 AS release

ARG TARGETARCH
ENV DORADO_VERSION=2.2.0

RUN test "${TARGETARCH}" = "amd64"

# SHA-256 of ONT's versioned Linux x64 release archive, verified over HTTPS.
ADD --checksum=sha256:c2376701a48bc6227a9cb5ad1d0b9ada60b5f3c8760a657066131fdf0d9545c4 \
    https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-linux-x64.tar.gz /tmp/dorado.tar.gz

RUN tar -xzf /tmp/dorado.tar.gz -C /opt \
    && mv "/opt/dorado-${DORADO_VERSION}-linux-x64" /opt/dorado \
    && rm /tmp/dorado.tar.gz

FROM nvidia/cuda:12.8.1-cudnn-runtime-ubuntu24.04@sha256:9175fa92f96de35a8cfb9493f0dfcf9435c7a597e9d95ad41d2cae382a95e3f9

ARG VCS_REF=unknown
ENV DORADO_VERSION=2.2.0 \
    PATH="/opt/dorado/bin:${PATH}"

LABEL org.opencontainers.image.title="Dorado 2.2.0" \
    org.opencontainers.image.version="2.2.0" \
    org.opencontainers.image.source="https://github.com/dbuckleychla/nf-core-oncoseq-leuk" \
    org.opencontainers.image.revision="${VCS_REF}"

COPY --from=release /opt/dorado /opt/dorado
WORKDIR /opt

# Keep the NVIDIA entrypoint shell-compatible for Nextflow task scripts.
CMD ["dorado"]
