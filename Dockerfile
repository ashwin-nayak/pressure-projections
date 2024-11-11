FROM julia:1.11.0

ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=${USER_UID}

# copy relevant project source files into container
ARG PROJECT_DIR=/pressure-projections
COPY ./ ${PROJECT_DIR}/

# Install dependencies
RUN julia --project=/pressure-projections -e 'using Pkg; Pkg.instantiate()'

# provide a non-root user for devcontainer
# the args for UID and GID will be overridden by VSCode
# Ref https://code.visualstudio.com/remote/advancedcontainers/add-nonroot-user
RUN addgroup --gid ${USER_GID} ${USERNAME} && \
    adduser --disabled-password --gecos '' --uid ${USER_UID} --gid ${USER_GID} ${USERNAME} && \
    usermod --uid ${USER_UID} --gid ${USER_GID} ${USERNAME}

USER ${USERNAME}
