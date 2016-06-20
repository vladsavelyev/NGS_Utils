import hashlib
from os import environ
import socket
import re
from collections import OrderedDict
from os.path import join, basename, splitext

from logger import info, critical, err, debug
from file_utils import file_exists, verify_file, file_transaction, adjust_path


class OrderedDefaultDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default

    def __reduce__(self):  # optional, for pickle support
        args = (self.default_factory,) if self.default_factory else ()
        return self.__class__, args, None, None, self.iteritems()


def _tryint(s):
    try:
        return int(s)
    except ValueError:
        return s


def _alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [_tryint(c) for c in re.split('([0-9]+)', s)]


def human_sorted(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=_alphanum_key)
    return l


def get_ext_tools_dirname(is_common_file=False):
    from sys import platform as _platform
    if not is_common_file and 'darwin' in _platform:
        return join('ext_tools', 'osx')
    else:
        return 'ext_tools'


def get_db_path(cnf, dbconf, dbname):
    db_path = cnf['genome'].get(dbname)
    if not db_path:
        db_path = dbconf.get('path')
        if not db_path:
            err('Please, provide a path to ' + dbname + ' in the "genomes" section in the system config. The config is: ' + str(cnf['genome']))
            return None
    return verify_file(db_path, is_critical=True)


def format_integer(name, value, unit=''):
    value = int(value)
    if value is not None:
        return '{name}: {value:,}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def format_decimal(name, value, unit=''):
    if value is not None:
        return '{name}: {value:.2f}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def mean(values):
    return float(sum(values)) / len(values) if len(values) > 0 else float('nan')


def median(values):
    values = sorted(values)

    if len(values) % 2 == 1:  # odd number of values
        return values[(len(values) - 1) / 2]
    else:  # even number of values - take the avg of central
        return (values[len(values) / 2] + values[len(values) / 2 - 1]) / 2


def get_numeric_value(string_value):
    """ parses string_value and returns only number-like part
    """
    num_chars = ['.', '+', '-']
    number = ''
    for c in string_value:
        if c.isdigit() or c in num_chars:
            number += c
    return number


def is_uk():
    hostname = socket.gethostname()
    return 'ukap' in hostname

def is_china():
    hostname = socket.gethostname()
    return 'cniclhpc' in hostname

def is_local():
    hostname = socket.gethostname()
    return 'local' in hostname or 'Home' in hostname or environ.get('PYTHONUNBUFFERED')

def is_us():
    hostname = socket.gethostname()
    return any(name in hostname for name in ['rask', 'chara', 'blue', 'green', 'espo', 'orr', 'usbod', 'bn0'])

def is_az():
    return is_us() or is_uk() or is_china() or is_sweden()

def is_cloud():
    hostname = socket.gethostname()
    return 'starcluster' in hostname

def is_sweden():
    hostname = socket.gethostname()
    return 'seml' in hostname


def md5(fpath):
    hash = hashlib.md5()
    with open(fpath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash.update(chunk)
    return hash.hexdigest()


def gray(text):
    return '<span class="gray">' + text + '</span>'
