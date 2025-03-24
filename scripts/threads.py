import concurrent.futures
from threading import Thread
import threading



def create_pool(n_threads):

    pool = concurrent.futures.ThreadPoolExecutor(max_workers=n_threads)
    return pool

def merge_pool(pool):
    pool.shutdown(wait=True)

def submit(pool, fun, *args, **kwargs):
    return pool.submit(fun, *args, **kwargs)


def create_thread(fun, *args, **kwargs):
    #print(ts, fun, kwargs)

    thread = Thread(target=thread_housing,  args=[fun] + list(args), kwargs=kwargs)
    return thread

def run_threads(list):
    for thread in list:
        thread.start()
        print("Launched:", thread.name)

def wait_threads():
    import time
    while threading.active_count() > 1:
        print("Alive threads:", threading.active_count(), threading.enumerate(), end ="\r")
        time.sleep(1)


def thread_housing(*args, **kwargs):
    fun = args[0]
    fun(*args[1:], **kwargs)
    thread = threading.current_thread()
    print("Finished:", thread.name)

