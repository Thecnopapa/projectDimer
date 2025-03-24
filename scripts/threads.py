import concurrent.futures




def create_pool(n_threads):

    pool = concurrent.futures.ThreadPoolExecutor(max_workers=n_threads)
    return pool

def merge_pool(pool):
    pool.shutdown(wait=True)

def submit(pool, fun, *args, **kwargs):
    return pool.submit(fun, *args, **kwargs)