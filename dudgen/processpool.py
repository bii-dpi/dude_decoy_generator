"""Easily manage either a pool or a FIFO queue of subprocesses.

Michael Mysinger 200709 Created
Michael Mysinger 200802 Added ProcessPool
"""

from collections import deque
import time

class ProcessQueue:
    """Manages a FIFO process queue of fixed size."""
    def __init__(self, size):
        self.size = size
        self.q = deque()

    def __iter__(self):
        """Support the iterator protocol."""
        return self

    def next(self):
        """Wait for old process and return its item."""
        try:
            process, item = self.q.popleft()
        except IndexError:
            raise StopIteration
        process.wait()
        return item

    def run(self, process, item=True):
        """Add process to queue and return old process item if full."""
        self.q.append((process, item))
        if len(self.q) >= self.size:
            return self.next()
        else:
            return None

class ProcessPool(ProcessQueue):
    """Manages a process pool where any running process can finish first."""
    def next(self):
        """Find a finished process and return its item."""
        delay =0.0005
        while True:
            # Find finished process if avaialble
            for i in xrange(self.size):
                try:
                    process, item = self.q[i]
                except IndexError:
                    if i == 0:
                        raise StopIteration
                    else:
                        break
                if process.poll() is not None:
                    del self.q[i]
                    return item
            # Or check again after waiting (exponential delay, max 1 sec)
            delay = min(delay*2, 1.0)
            time.sleep(delay)
