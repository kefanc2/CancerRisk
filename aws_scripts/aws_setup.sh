wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh

bash Anaconda3-2020.02-Linux-x86_64.sh

export PATH=~/anaconda3/bin:$PATH

echo 'export PATH=~/anaconda3/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

conda create -y -n py2 python=2.7

sudo yum install -y gcc-c++

sudo yum install -y gsl gsl-devel

sudo yum install tmux

cd phylowgs

g++ -o mh.o -O3 mh.cpp  util.cpp `gsl-config --cflags --libs`




