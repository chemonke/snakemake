# Stage 1: Build Surge
FROM ubuntu:22.04 AS build

# Install required build tools and dependencies
RUN apt-get update && \
    apt-get -y install curl gcc make zlib1g-dev && \
    apt-get -y clean && \
    rm -rf \
      /var/lib/apt/lists/* \
      /usr/share/doc \
      /usr/share/doc-base \
      /usr/share/man \
      /usr/share/locale \
      /usr/share/zoneinfo

# Set working directory
WORKDIR /root

# Download and build nauty
RUN curl -L -o nauty27r3.tar.gz https://users.cecs.anu.edu.au/~bdm/nauty/nauty27r3.tar.gz && \
    tar xzvf nauty27r3.tar.gz && \
    cd nauty27r3 && \
    ./configure && make

# Set environment variable for nauty
ENV NAUTY_HOME=/root/nauty27r3

# Copy Surge source files
COPY src/surge.c $NAUTY_HOME
COPY src/Makefile /root

# Build Surge
WORKDIR $NAUTY_HOME
RUN ln -s /root/nauty27r3 /root/nauty && \
    make -f ../Makefile clean && make -f ../Makefile surge

# Stage 2: Runtime environment
FROM ubuntu:22.04

# Install runtime dependencies
RUN apt-get update && \
    apt-get -y install curl time gnupg zlib1g && \
    apt-get -y clean && \
    rm -rf \
      /var/lib/apt/lists/* \
      /usr/share/doc \
      /usr/share/doc-base \
      /usr/share/man \
      /usr/share/locale \
      /usr/share/zoneinfo

# Copy Surge binary from the build stage
COPY --from=build /root/nauty27r3/surge /usr/bin
