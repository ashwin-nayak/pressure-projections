FROM julia:1.11.0

ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=${USER_UID}
ENV JULIA_DEPOT_PATH="/tmp"

RUN addgroup --gid ${USER_GID} ${USERNAME} && \
    adduser --disabled-password --gecos '' --uid ${USER_UID} --gid ${USER_GID} ${USERNAME} && \
    usermod --uid ${USER_UID} --gid ${USER_GID} ${USERNAME}

USER ${USERNAME}
