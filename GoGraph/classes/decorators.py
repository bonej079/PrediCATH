import functools
import logging

import logzero
from logzero import logger
from logzero import setup_logger


# "/tmp/phdlogs/main-logfile.log"
def create_logger(path: str = None, loglevel=logging.INFO):
    """
    Creates a logging object and returns it
    """
    logzero.loglevel(loglevel)

    # create the logging file handler
    if path is not None:
        import os
        if os.path.isdir(path):
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)
            path = path.rstrip(os.path.sep)
            path += os.path.sep + "main-logfile.log"

        if not os.path.exists(path):
            directories = os.path.dirname(path)
            if not os.path.exists(directories):
                os.makedirs(directories)
            os.mknod(path)

        logzero.logfile(path, maxBytes=100000000, backupCount=100)  # 100 MB * 10 = 1GB

    fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(fmt)
    logzero.formatter(formatter)

    # add handler to logger object
    return logger


# "/tmp/phdlogs/main-logfile.log"
def create_custom_logger(path: str = None, name="default_logger", loglevel=logging.INFO):
    """
    Creates a logging object and returns it
    """
    fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(fmt)
    logger = setup_logger(name=name, logfile=path, level=loglevel, maxBytes=100000000, backupCount=100, formatter=formatter)

    # create the logging file handler
    if path is not None:
        import os
        if os.path.isdir(path):
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)
            path = path.rstrip(os.path.sep)
            path += os.path.sep + "main-logfile.log"

        if not os.path.exists(path):
            directories = os.path.dirname(path)
            if not os.path.exists(directories):
                os.makedirs(directories)
            os.mknod(path)

    # add handler to logger object
    return logger

excLogger = create_logger()

def exception(excLogger):
    """
    A decorator that wraps the passed in function and logs
    exceptions should one occur

    @param logger: The logging object
    """

    def decorator(function):
        @functools.wraps(function)
        def wrapper(*args, **kwargs):
            try:
                return function(*args, **kwargs)
            except:
                # log the exception
                err = "There was an exception in  "
                err += function.__name__
                excLogger.exception(err)

                # re-raise the exception
                raise
        return wrapper
    return decorator
