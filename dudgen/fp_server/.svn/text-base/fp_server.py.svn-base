#!/usr/arch/bin/python -u
# Fingerprint-generating XML-RPC server
#
# Michael Keiser   200607 Created
# Michael Mysinger 200608 Modified to use fastdl extension module
# Michael Mysinger 200702 Modified for remote startup and shutdown via ssh
# Michael Mysinger 200808 Switch to static Forking Mixin server

import sys
import socket
from SimpleXMLRPCServer import SimpleXMLRPCServer
from SocketServer import ForkingMixIn
import fastdl

class ForkingXMLRPCServer(ForkingMixIn, SimpleXMLRPCServer):
        pass

QUIT = False

# smiles = [(smi, key), ...]
# return = [(asciifp, key), ...]
def get_fingerprints(smiles):
    return fastdl.fingerprints(smiles)

# main
if __name__ == '__main__':

    if len(sys.argv) == 2:
        port = int(sys.argv[1])
    else:
        port = int(raw_input('port number: '))

    print 'Starting fingerprint server on %d.' % port
    socket.setdefaulttimeout(60)
    server = ForkingXMLRPCServer(('', port), allow_none=True)
    server.register_function(get_fingerprints)

    server.serve_forever()
