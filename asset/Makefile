SHELL=/bin/bash -e -x -o pipefail
all:
	
	sudo apt-get update
	sudo apt-get --yes install git
	sudo git clone https://github.com/chrchang/eigensoft.git
	cp eigensoft/bin/* /usr/bin
	
	# install GNU Scientific library
	sudo apt-get --yes install gsl-bin libgsl0-dbg libgsl0-dev libgsl0ldbl
	# SQLite
	sudo apt-get --yes install sqlite3 libsqlite3-dev
	
	sudo apt-get --yes install pigz tabix openjdk-8-jre-headless samtools picard-tools parallel vcftools 
	
