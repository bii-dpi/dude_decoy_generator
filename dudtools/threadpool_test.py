import threadpool
import time

def print_thread_name(self):
    print '\n' + self.getName()
    time.sleep(0.1)

def nested_loop(size):
    return sum((sum((x*y for x in xrange(size))) for y in xrange(size)))

def print_result(req, result):
    print 'The size was', req.args[0], 'and the result was', result

def exc_test():
    print 'Throwing ValueError exception in thread...'
    raise ValueError('A value error occurred in this thread.')

def exc_print(req, error):
    print 'Processing exception callback in main...'
    print 'Type: ', error[0]
    print 'Value: ', error[1]

tp =  threadpool.ThreadPool(5)
preq = threadpool.Request(print_thread_name, sendself=True)
for i in xrange(50):
    tp.put(preq)

time.sleep(3)
tp.poll()
tp.put(threadpool.Request(exc_test, exc_callback=exc_print))
time.sleep(5)
tp.get()
tp.shutdown()
del tp

tp =  threadpool.ThreadPool(2, qsize=10)
for i in xrange(18):
    tp.put(threadpool.Request(nested_loop, (1500,), callback=print_result))

tp.poll()
tp.poll()
tp.get()
tp.get()
tp.get()
time.sleep(4)
tp.poll()
tp.wait()
tp.shutdown()
del tp
