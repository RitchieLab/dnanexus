#!/usr/bin/env python

import sys
import os
#import multiprocessing
import dxpy
import io


class randf(io.RawIOBase):

    def read(self, n):
        return os.urandom(n)

    def fileno(self):
        return 0
    
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
    
    # create a random named file
    
    if len(sys.argv) > 1:
        given_rpath=sys.argv[2]
        if given_rpath.startswith("/"):
            remote_path=given_rpath
        else:
            remote_path = remote_path + "/" + given_rpath

    f = randf()
    
    dxpy.bindings.dxfile_functions.upload_local_file(file=f, folder=remote_path, parents=True, project=curr_proj)

