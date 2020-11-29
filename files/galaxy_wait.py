#!/usr/bin/env python
'''Utility to do a blocking sleep until a Galaxy instance is responsive.
This is useful in docker images, in RUN steps, where one needs to wait
for a currently starting Galaxy to be alive, before API requests can be
made successfully.

The script functions by making repeated requests to
``http(s)://fqdn/api/version``, an API which requires no authentication
to access.'''
from argparse import ArgumentParser
import requests
import sys
import time






def _parser():
    '''Constructs the parser object'''
    parser = ArgumentParser( usage="usage: %(prog)s <options>",
                            description="Script to sleep and wait for Galaxy to be alive.")
    parser.add_argument("--timeout",
                        default=60, type=int,
                        help="Galaxy startup timeout in seconds. The default value is 60")
    parser.add_argument("--galaxy",
                        default="http://127.0.0.1:8080", 
                        help="Galaxy URL")
    return parser


def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = _parser()
    return parser.parse_args()


def galaxy_wait(galaxy_url, timeout=600, verbose=False):
    count = 0
    while True:
        try:
            result = requests.get(galaxy_url + '/api/version')
            try:
                result = result.json()
                if verbose:
                    sys.stdout.write("Galaxy Version: %s\n" % result['version_major'])
                break
            except ValueError:
                if verbose:
                    sys.stdout.write("[%02d] No valid json returned... %s\n" % (count, result.__str__()))
                    sys.stdout.flush()
        except requests.exceptions.ConnectionError as e:
            if verbose:
                sys.stdout.write("[%02d] Galaxy not up yet... %s\n" % (count, e[0:100]))
                sys.stdout.flush()
        count += 1

        # If we cannot talk to galaxy and are over the timeout
        if timeout != 0 and count > timeout:
            sys.stderr.write("Failed to contact Galaxy\n")
            sys.exit(1)

        time.sleep(1)


def main():
    """
    Main function
    """
    options = _parse_cli_options()

    galaxy_wait(options.galaxy,options.timeout)

    sys.exit(0)


if __name__ == "__main__":
    main()
