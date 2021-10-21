"""General thread pool using input and output queues and allowing callbacks 

Request class
Wrap a function defining a unit of work for the thread pool.

GeneralThread class
Background worker thread servicing input requests to generate output.

ThreadPool class
Thread pool distributes queue requests and processes the results.

make_requests helper function
Generate a list of Requests for target differing only in the arguments.

Michael Mysinger 200608 Created
Michael Mysinger 200803 Switch to new-style classes
"""

import sys
import Queue
import socket
import threading
import xmlrpclib
from collections import deque

class Request(object):
    """Wrap a function defining a unit of work for the thread pool.

    target is the function to wrap.
    args is a tuple containing the arguments to the function.
    kwargs is a dictionary containing the keyword arguments.
    sendself is a boolean to request that thread object be sent to target.
    callback is called in the main thread when this request is processed.
    exc_callback is called in the main thread if the their was an exception.
    g is used to hold a user object for access to globals in a callback chain
    
    A Request object passed to the thread pool is included with the result.
    Thus all arguments are available when you process the data. In addition,
    you can add needed per-request data to the Request object so you may
    easily access it later.
    """
    def __init__(self, target, args=None, kwargs=None, sendself=None, 
                 callback=None, exc_callback=None, g=None):
        self.exception = False
        self.target = target
        if args is None:
            args = []
        self.args = args
        if kwargs is None:
            kwargs = {}
        self.kwargs = kwargs
        self.callback = callback
        self.exc_callback = exc_callback
        self.sendself = sendself
        self.g = g

class GeneralThread(threading.Thread):
    """Background worker thread servicing input requests to generate output.

    num is the number assigned to this thread.
    inqueue is the thread pool input queue.
    outqueue is the thread pool output queue.
    init is called as the thread is initialized.
    cleanup is called as the thread exits.

    Together init and cleanup (with Request.sendself) allow storage of
    (and access to) per-thread data.
    """
    def __init__(self, num, inqueue, outqueue, init=None, cleanup=None):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.inqueue = inqueue
        self.outqueue = outqueue
        self.num = num
        self.cleanup = cleanup
        if init:
            assert isinstance(init, Request)
            if init.sendself:
                init.target(self, *init.args, **init.kwargs)
            else:
                init.target(*init.args, **init.kwargs)
        self.start()

    def run(self):
        """Waits for input requests and generates output results."""
        while True:
            req = self.inqueue.get()
            if req is None:
                cleanup = self.cleanup
                if cleanup:
                    assert isinstance(cleanup, Request)
                    if cleanup.sendself:
                        cleanup.target(self, *cleanup.args, **cleanup.kwargs)
                    else:
                        cleanup.target(*cleanup.args, **cleanup.kwargs)
                    
                break
            else:
                try:
                    if req.sendself:
                        result = req.target(self, *req.args, **req.kwargs)
                    else:
                        result = req.target(*req.args, **req.kwargs)
                    self.outqueue.put((req, result))
                except:
                    req.exception = True
                    self.outqueue.put((req, sys.exc_info()))

class ThreadPool(object):
    """Thread pool to distribute queue requests and process the results.

    num_threads is the number of threads in the pool.
    init is called during initialization of each member thread.
    cleanup is called by each member thread during thread pool shutdown.
    qsize lets you limit the size of the input queue.

    Together init and cleanup (with Request.sendself) allow storage of
    (and access to) per-thread data.
    """
    def __init__(self, num_threads, init=None, cleanup=None, qsize=0):
        self.pool = []
        self.inqueue = Queue.Queue(qsize)
        self.outqueue = Queue.Queue()
        self.qsize = 0 
        for i in xrange(num_threads):
            self.pool.append(
                GeneralThread(i, self.inqueue, self.outqueue,
                              init=init, cleanup=cleanup))

    def put(self, req, block=True):
        """Put new request into input queue (blocking if full by default).

        If you change block to False, be sure to catch Queue.Full.
        """
        assert isinstance(req, Request)
        self.inqueue.put(req, block)
        self.qsize += 1

    def get(self, block=True):
        """Process next result in the queue (blocking for it by default).

        If you change block to False, be sure to catch Queue.Empty (see poll).
        """
        req, result = self.outqueue.get(block)
        self.qsize -= 1
        if req.exception:
            if req.exc_callback:
                req.exc_callback(req, result)
            else:
                exc_type, exc_value, exc_traceback = result
                raise exc_type, exc_value, exc_traceback
        elif req.callback:
            req.callback(req, result)
        else:
            return req, result

    def poll(self):
        """ Process all output results currently waiting in the queue.""" 
        results = []
        while True:
            try:
                result = self.get(False)
                if result is not None:
                    results.append(result)
            except Queue.Empty:
                break
        return results or None

    def wait(self):
        """ Wait for all results to finish, blocking until they arrive.""" 
        results = []
        for i in xrange(self.qsize):
            result = self.get()
            if result is not None:
                results.append(result)
        return results or None

    def shutdown(self):
        """ Shutdown thread pool by stopping all my threads.""" 
        for i in self.pool:
            self.inqueue.put(None)
        for i in self.pool:
            i.join()

def make_requests(target, args_list, kwargs_list=None,
                 callback=None, exc_callback=None, sendself=None):
    """Generate a list of Requests for target differing only in arguments."""
    args_len = len(args_list)
    if kwargs_list and args_len != len(kwargs_list):
        raise IndexError(
            'If given, kwargs_list must be same length as args_list!')

    requests = []
    for i in xrange(args_len):
        if kwargs_list:
            requests.append(Request(target, args_list[i], kwargs_list[i],
                                    callback, exc_callback, sendself))
        else:
            requests.append(Request(target, args_list[i], None,
                                    callback, exc_callback, sendself))

    return requests

class OverflowTP(ThreadPool):
    """Class to help manage out-of-order threadpool results. 

    Although, things go into a threadpool input queue in order, they often
    come out in a different order due to parallel processing and variable
    task size. This threadpool extension helps manage this unpredicitable
    order. Early results can be returned to the overflow queue of the class
    using the putback method. When the reset method is called, previously
    overflowed items are returned to the head of the output results, and
    are preferentially returned by the get method.

    A common use is to reorder results back to the normal nested loop
    order that was originally used to insert items into the queue.
    """
    def __init__(self, num_threads, init=None, cleanup=None, qsize=0):
        ThreadPool.__init__(self, num_threads, init, cleanup, qsize)
        self.overflows = deque()
        self.overflowed = deque()

    def get(self, block=True):
        """Get results from the previous overflow queue preferentially."""
        if len(self.overflowed) > 0:
            return self.overflowed.popleft()
        else:
            return ThreadPool.get(self, block)

    def putback(self, req, res):
        """Put back early results into the overflow queue."""
        self.overflows.append((req, res)) 

    def reset(self, extras=None):
        """Reset the overflow queue to empty and get previous overflows."""
        self.overflowed = self.overflows
        self.overflows = deque()
        if extras:
            self.overflowed.extend(extras)

def init_client(self, host, port):
    """Establish local xmlrpc proxy."""
    self.client = xmlrpclib.ServerProxy('http://%s:%d' % (host, port), 
                                        allow_none=True)

class RemotePool(ThreadPool):
    """Local threadpool connected a series of remote XML-RPC servers. 

    Uses  XML-RPC to transmit and receive data from remote servers, 
    and a local ThreadPool to send work requests to them. Can 
    optionally use ssh to dynamically start the remote server.
    """
    def __init__(self, num_threads, host, port, script, qsize=0, spawn=False):
        self.host = host
        self.port = port
        self.script = script
        if spawn:
            import pexpect
            ssh = ['ssh', host, script, "0.0.0.0", str(port)]
            self.server = pexpect.spawn(' '.join(ssh))
        init = Request(init_client, (host, port), sendself=True)
        ThreadPool.__init__(self, num_threads, init=init, qsize=qsize)

