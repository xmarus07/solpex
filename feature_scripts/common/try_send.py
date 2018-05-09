import requests as r
import time
import sys


def try_send(fun, args, default="", max_tries=5, wait_fun=lambda x: x * x + 3):
    """
    Try execution of function fun several times until it succeed or max tries
    is reached. If max tries is reached default value will be returned.
    """
    tries = 0
    while True:
        try:
            return fun(args)
        except r.ConnectionError:
            print("Connection Error.", file=sys.stderr, flush=True)
            time.sleep(wait_fun(tries))
        except Exception as e:
            print(e.args, file=sys.stderr)
            return default
        finally:
            tries = tries + 1
            if tries == max_tries:
                return default
