#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get update
apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  build-essential \
  cmake \
  automake \
  gcc \
  gcc-10 \
  g++-10 \
  gdb \
  git \
  python3 \
  python3-pip \
  python3-setuptools \
  python3-dev \
  wget \
  zlib1g-dev \
  libbz2-dev \
  libgomp1 \
  liblzma-dev \
  libhts-dev \
  tabix \
  curl \
  libvcflib-tools \
  libcurl4-gnutls-dev \
  libssl-dev \
  samtools

update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 --slave /usr/bin/g++ g++ /usr/bin/g++-10 --slave /usr/bin/gcov gcov /usr/bin/gcov-10
python3 -m pip install tox

if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root

#_____________________ enaBrowserTools ____________________#
cd $install_root
git clone https://github.com/enasequence/enaBrowserTools.git


#_________________________ NGmerge __________________________#
cd $install_root
git clone https://github.com/harvardinformatics/NGmerge.git NGmerge-git
cd NGmerge-git
git checkout 224fc6a0066024e05965d101d998704815cb4c41
make
cd ..
mv NGmerge-git/NGmerge .
rm -rf NGmerge-git


#_________________________ bcftools _________________________#
cd $install_root
wget -q https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar xf bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2/
make
cd ..
mv bcftools-1.10.2/bcftools .
rm -rf bcftools-1.10.2

#________________________ vt ________________________________#
cd $install_root
git clone https://github.com/atks/vt.git vt-git
cd vt-git
git checkout 2187ff6347086e38f71bd9f8ca622cd7dcfbb40c
make
cd ..
mv vt-git/vt .
rm -rf vt-git

#________________________ minimap2 __________________________#
cd $install_root
git clone https://github.com/lh3/minimap2.git minimap2_git
cd minimap2_git
git checkout b0b199f5039e8da5e8d1a9a7ae130580fd33fe1f
arch_is_arm=$(dpkg --print-architecture | grep '^arm' | wc -l)
if [[ $arch_is_arm -gt 0 ]]
then
    make arm_neon=1 aarch64=1
else
    make
fi
cd ..
mv minimap2_git/minimap2 .
rm -rf minimap2_git

#________________________ racon _____________________________#
cd $install_root
git clone --recursive https://github.com/lbcb-sci/racon.git racon-git
cd racon-git
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make CC=gcc-10 CPP=g++-10 CXX=g++-10 LD=g++-10
cd ../../
mv racon-git/build/bin/racon .
rm -rf racon-git

#________________________ mummer ____________________________#
cd $install_root
wget -q https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -xvf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure LDFLAGS=-static
make
make install
ldconfig
cd ..
rm -rf mummer-4.0.0rc1

#________________________ cylon _____________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/cylon.git
cd cylon
git checkout 8b61712aaff9a674f497b6569312cd84bfc446b6
pip3 install .
cd ..
rm -rf cylon

#______________________ ReadItAndKeep _______________________#
cd $install_root
git clone https://github.com/GlobalPathogenAnalysisService/read-it-and-keep.git
cd read-it-and-keep
git checkout 61ae15be1e515c960b0135eae7dd59568a9de30d
cd src
make
cd $install_root
mv read-it-and-keep/src/readItAndKeep .
rm -rf read-it-and-keep

#_________________________ mafft __________________________#
# Can't apt get mafft because the version is too old and doesn't work
# with how we call from python (avoid using files). Install from source.
# See https://mafft.cbrc.jp/alignment/software/installation_without_root.html
cd $install_root
wget https://mafft.cbrc.jp/alignment/software/mafft-7.525-without-extensions-src.tgz
tar xf mafft-7.525-without-extensions-src.tgz
cd mafft-7.525-without-extensions/core
sed  -i "s~PREFIX = /usr/local~PREFIX = $install_root/mafft_install~" Makefile
make
make install
cd $install_root
cp -s mafft_install/bin/mafft .
rm -rf mafft-7.525-without-extensions*

#________________________ varifier __________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/varifier.git
cd varifier
git checkout 940d4503671f5dd09ce51116be6e77c85d50ada5
pip3 install .
cd ..
rm -rf varifier


apt-get remove -y \
  build-essential \
  cmake \
  automake \
  gcc \
  gcc-10 \
  g++-10 \
  gdb \
  git \
  wget \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libhts-dev \
  tabix \
  curl \
  libvcflib-tools \
  libcurl4-gnutls-dev \
  libssl-dev
