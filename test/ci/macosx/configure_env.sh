export SPACK_DISABLE_LOCAL_CONFIG=true
export SPACK_USER_CACHE_PATH=".spack-cache"

# spack environment
git clone -b v0.20.2 https://github.com/spack/spack
source spack/share/spack/setup-env.sh

mkdir spack/var/spack/environments
cp -r ./test/ci/macosx spack/var/spack/environments

spack env activate macosx
spack install
