export SPACK_DISABLE_LOCAL_CONFIG=true
export SPACK_USER_CACHE_PATH=".spack-cache"

SPACK_RELEASE=v0.20.2

# spack environment
if [ ! -d spack ]; then
  git clone -b ${SPACK_RELEASE} https://github.com/spack/spack
else
  cd spack
  git checkout ${SPACK_RELEASE}
  cd -
fi
source spack/share/spack/setup-env.sh

mkdir spack/var/spack/environments
cp -r ./test/ci/macosx spack/var/spack/environments

spack env activate macosx
spack install
