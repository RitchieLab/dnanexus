SHELL=/bin/bash -e -x -o pipefail
all:
	
	sudo apt-get update
	sudo apt-get --yes install git
	sudo git clone https://github.com/chrchang/eigensoft.git
	sudo cp eigensoft/bin/* /usr/bin
	
	sudo curl -L https://cpanmin.us/ -o cpanm
	sudo chmod +x cpanm && sudo mv cpanm /usr/bin
	
	sudo pip install apsw	
	
	sudo apt-get --yes install liblzma-dev pigz tabix openjdk-8-jre-headless parallel vcftools	
	sudo apt-get --yes install libtest-mockmodule-perl libcpan-meta-yaml-perl libcpan-meta-perl libdbi-perl r-recommended	
	# install GNU Scientific library
	sudo apt-get --yes install gsl-bin libgsl0-dbg libgsl0-dev libgsl0ldbl
	# SQLite
	sudo wget https://www.sqlite.org/2018/sqlite-autoconf-3240000.tar.gz 
	sudo tar -zxvf sqlite-autoconf-3240000.tar.gz 
	pushd sqlite-autoconf-3240000 && ./configure --enable-dynamic-extensions && sudo make && sudo make install && popd 
	sudo wget -c "https://sqlite.org/contrib/download/extension-functions.c/download/extension-functions.c?get=25" -O extension-functions.c
	sudo gcc -fPIC -lm -shared extension-functions.c -o libsqlitefunctions.so && mv libsqlitefunctions.so /usr/local/lib
	# Picard tools
	sudo wget https://github.com/broadinstitute/picard/releases/download/2.18.9/picard.jar && mv picard.jar /
	# HTSlib
	sudo wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2 && tar -xjf htslib-1.8.tar.bz2
	pushd htslib-1.8 && sudo ./configure && sudo make && sudo make install && popd
	# bcftools
	sudo wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2 && tar -xjf bcftools-1.8.tar.bz2
	pushd bcftools-1.8 && sudo ./configure --prefix=/usr/local/ && sudo make -s && sudo make install && popd
	# samtools	
	sudo wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && tar xfj samtools-1.7.tar.bz2
	pushd samtools-1.7 && sudo ./configure && sudo make && sudo make install && popd
