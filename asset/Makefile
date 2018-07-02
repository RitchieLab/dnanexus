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
	sudo apt-get --yes install sqlite3 libsqlite3-dev
	
	sudo apt-get --yes install pigz tabix openjdk-8-jre-headless samtools picard-tools parallel vcftools 
	sudo apt-get --yes install libtest-mockmodule-perl libcpan-meta-yaml-perl libcpan-meta-perl libdbi-perl r-recommended	
