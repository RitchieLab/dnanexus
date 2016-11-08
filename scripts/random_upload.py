#!/usr/bin/env python

import sys
import os
#import multiprocessing
import dxpy
import io
import socket


class randf(io.RawIOBase):

    def read(self, n):
        return os.urandom(n)

    def fileno(self):
        return 0
    
if __name__ == "__main__":

    upload_path="."
    dx_path=""
    remote_path="/"
    curr_proj = "project-Bz06Zv80pjKvypFvX4kqp1Vk"
    API_token={"auth_token_type" : "Bearer", "auth_token" : "Y1WoDbANKOKDF0ZfPV0wt9cQaIotIOWC"}
    # try to find the remote path
    # by default use "/From_<hostname>"
    remote_path = "/From_" + socket.gethostname()
    
    # create a random named file
    
    if len(sys.argv) > 1:
        given_rpath=sys.argv[1]
        if given_rpath.startswith("/"):
            remote_path=given_rpath
        else:
            remote_path = remote_path + "/" + given_rpath
  
    # set the security context
    dxpy.set_security_context(API_token)
    
    # upload to infinity and beyond!
    f = randf()
    dxpy.bindings.dxfile_functions.upload_local_file(file=f, folder=remote_path, parents=True, project=curr_proj)

