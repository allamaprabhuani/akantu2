FROM manjarolinux/base:latest as base

# squashing the whole base image into one layer
FROM scratch AS release
COPY --from=base / /

ARG TARGETPLATFORM

ARG CACHEBUST=1
ENV LANG=en_US.UTF-8

ENV PATH="/usr/bin:${PATH}"

RUN uname -m && \
    pacman-key --init && \
    pacman-mirrors -f 5

RUN [[ "${TARGETPLATFORM}" == "linux/amd64" ]] || exit 0 && \
    pacman -Syy --noconfirm --needed archlinux-keyring manjaro-keyring && \
    pacman-key --populate archlinux manjaro

RUN [[ "${TARGETPLATFORM}" == "linux/arm64" ]] || exit 0 && \
    pacman -Syy --noconfirm --needed archlinuxarm-keyring manjaro-arm-keyring && \
    pacman-key --populate archlinuxarm manjaro-arm

RUN pacman -S --noconfirm --needed \
  shadow \
  git \
  git-lfs \
  cmake \
  libseccomp \
  autoconf \ 
  automake \
  binutils \
  bison  \
  fakeroot \
  file \
  findutils \
  flex \
  gawk \
  gcc \
  gettext \
  grep \
  groff \
  gzip \
  libtool \
  m4 \
  make \
  pacman \
  patch \
  pkgconf \
  sed  \
  sudo \
  texinfo \
  lsb-release \
  manjaro-release \
  boost \
  which && \
  # docker context give real space limits
  sed -i -e 's~CheckSpace.*~#CheckSpace~g' '/etc/pacman.conf' && \
  pacman -Syyu --noconfirm --needed

RUN ls /etc/*-release && cat /etc/*-release

# user 'builder' can be used as the running user for applications prohibiting root usage (pacman)
RUN id -u builder &>/dev/null || (useradd -d /builder -m builder && \
    echo "builder ALL=(ALL) NOPASSWD: ALL" >>/etc/sudoers)

RUN pacman -Syu --noconfirm rust

USER builder

RUN cd && git clone https://aur.archlinux.org/paru.git && \
  cd paru && makepkg -s && \
  sudo pacman -U --noconfirm paru*.zst && \
  cd /builder && rm -rf paru

RUN paru -S --noconfirm scotch scalapack parmetis metis

RUN cd && curl -O http://mumps.enseeiht.fr/MUMPS_5.4.1.tar.gz && \
  tar xvzf MUMPS_5.4.1.tar.gz

COPY Makefile.inc /builder/MUMPS_5.4.1

RUN cd /builder/MUMPS_5.4.1 && make

CMD ["/usr/bin/bash"]
