#!/usr/bin/env python

import sys
import os
import multiprocessing
import dxpy


class Uploader(object):
	def __init__(self, basedir, proj, rempath=""):
		self.basedir = basedir
		self.rempath = rempath
		self.proj = proj
	
	def __call__(self, fullf):
		folder_str=self.rempath + "/" + '/'.join((fullf[len(self.basedir):] if fullf.startswith(self.basedir) else fullf).split(os.path.sep)[:-1])
		folder_str=("/" + folder_str if not folder_str.startswith("/") else folder_str)
	
		print >> sys.stderr, "Uploading " + fullf + " to " + folder_str
		dxpy.bindings.dxfile_functions.upload_local_file(filename=fullf, folder=folder_str, parents=True, project=self.proj)
		print >> sys.stderr, "Finished uploading "+ fullf

if __name__ == "__main__":

	upload_path="."
	dx_path=""
	remote_path="/"
	curr_proj = dxpy.PROJECT_CONTEXT_ID
	# try to find the remote path
	try:
		with open(os.path.join(os.path.expanduser("~"), ".dnanexus_config", "DX_CLI_WD")) as f:
			remote_path=f.read()
	except:
		pass
	
	if len(sys.argv) > 1:
		upload_path=sys.argv[1]
	if len(sys.argv) > 2:
		given_rpath=sys.argv[2]
		if given_rpath.startswith("/"):
			remote_path=given_rpath
		else:
			remote_path = remote_path + "/" + given_rpath
		
	pathlist = []
	for r,d,fs in os.walk(upload_path):
		for f in fs:
			pathlist.append(os.path.join(r,f))
			
	
	p = multiprocessing.Pool(multiprocessing.cpu_count())
	p.map(Uploader(upload_path, curr_proj, remote_path), pathlist)
