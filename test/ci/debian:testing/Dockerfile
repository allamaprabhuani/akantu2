FROM debian:testing-slim

# for release generations
RUN apt -qq update && apt -qq -y install \
    cargo curl git \
    && rm -rf /var/lib/apt/lists/*

RUN cargo install markdown-extract \
    && mv ~/.cargo/bin/markdown-extract /usr/bin \
    && rm -rf ~/.cargo

RUN curl -fSsL -o /usr/bin/release-cli \
    https://gitlab.com/gitlab-org/release-cli/-/releases/v0.16.0/downloads/bin/release-cli-linux-amd64 \
    && chmod +x /usr/bin/release-cli
