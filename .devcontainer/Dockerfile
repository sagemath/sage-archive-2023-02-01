ARG VARIANT="jammy"
FROM ghcr.io/sagemath/sage/sage-docker-ubuntu-${VARIANT}-standard-with-system-packages:dev

# Install additional packages needed for devcontainer support in VS code
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends gpgconf openssh-client
