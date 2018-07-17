all:
	sudo apt-get install -y libtool
	# First change out boost
	sudo apt-get purge -y --auto-remove libboost-dev
	wget http://downloads.sourceforge.net/project/boost/boost/1.48.0/boost_1_48_0.tar.bz2
	tar jxvf boost_1_48_0.tar.bz2
	# Patch it to get threads to work
	cd boost_1_48_0/ && wget https://svn.boost.org/trac10/raw-attachment/ticket/6165/libstdcpp3.hpp.patch && patch -p0 < libstdcpp3.hpp.patch
	sed -i 's/TIME_UTC/TIME_UTC_/g' boost_1_48_0/boost/thread/xtime.hpp
	sed -i 's/TIME_UTC/TIME_UTC_/g' boost_1_48_0/libs/thread/src/pthread/timeconv.inl
	# And compile boost
	cd boost_1_48_0/ && ./bootstrap.sh --prefix=/usr/local --with-libraries=all --libdir=/usr/local/lib --includedir=/usr/local/include && sudo ./b2 install
	cd boost_1_48_0/ && ./bootstrap.sh --prefix=/usr/local --with-libraries=thread --libdir=/usr/local/lib --includedir=/usr/local/include && sudo ./b2 install
	sudo /bin/bash -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/boost.conf'
	sudo ldconfig

	sudo apt-get install -y libgsl0-dev

	wget https://dl.dnanex.us/F/D/77pyvj79g4g9g2BQjj9bZ392V9Y47b5PY6f5v3kV/plato-2.1.0-beta3.tar.gz
	tar -zxvf plato-2.1.0-beta3.tar.gz --no-same-owner
	cd plato-2.1.0-beta3; autoreconf -i; mkdir -p build; cd build; ../configure --with-boost-libdir=/usr/local/lib; make && sudo make install

