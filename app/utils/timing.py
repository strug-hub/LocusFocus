from time import perf_counter
from flask import current_app as app


def timeit(func):
    """
    Log to flask logger the time it takes to run a function.
    Timing is logged at debug level.
    """
    def wrapper(*args, **kwargs):
        t1 = perf_counter()
        result = func(*args, **kwargs)
        t2 = perf_counter()
        app.logger.debug(f"{func.__name__} took {t2 - t1:.5f} seconds")
        return result

    return wrapper
