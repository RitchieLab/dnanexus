#!/usr/bin/env python

import sys
import os
import multiprocessing
import dxpy


class Uploader(object):
	def __init__(self, basedir):
		self.basedir = basedir
	
	def __call__(self, fullf):
		folder_str='/'.join((fullf[len(self.basedir):] if fullf.startswith(self.basedir) else fullf).split(os.path.sep)[:-1])
		folder_str=("/" + folder_str if not folder_str.startswith("/") else folder_str)
	
		dxpy.bindings.dxfile_functions.upload_local_file(filename=fullf, folder=folder_str, parents=True)

if __name__ == "__main__":

	upload_path="."
	dx_path=""
	if len(sys.argv) > 1:
		upload_path=sys.argv[1]
		
	pathlist = []
	for r,d,fs in os.walk(upload_path):
		for f in fs:
			pathlist.append(os.path.join(r,f))
			
	
	p = multiprocessing.Pool(multiprocessing.cpu_count())
	p.map(Uploader(upload_path), pathlist)
