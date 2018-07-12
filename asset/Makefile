SHELL=/bin/bash -e -x -o pipefail
all:
	
	sudo apt-get update
	sudo apt-get --yes install git
	sudo git clone https://github.com/chrchang/eigensoft.git
	sudo cp eigensoft/bin/* /usr/bin
	
	sudo curl -L https://cpanmin.us/ -o cpanm
	sudo chmod +x cpanm && sudo mv cpanm /usr/bin
	
	sudo pip install apsw	
	
	# install GNU Scientific library
	sudo apt-get --yes install gsl-bin libgsl0-dbg libgsl0-dev libgsl0ldbl
	# SQLite
	sudo apt-get --yes install libsqlite3-0
	# Picard tools
	sudo wget https://github.com/broadinstitute/picard/releases/download/2.18.9/picard.jar && mv picard.jar /
	
	sudo apt-get --yes install pigz tabix openjdk-8-jre-headless parallel vcftools 
	sudo apt-get --yes install liblzma-dev; wget  "https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2"  && tar xfj samtools-1.7.tar.bz2  && (cd samtools-1.7 && make)
	sudo apt-get --yes install libtest-mockmodule-perl libcpan-meta-yaml-perl libcpan-meta-perl libdbi-perl r-recommended	
