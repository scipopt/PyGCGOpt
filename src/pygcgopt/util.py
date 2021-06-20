import sys

if sys.version_info >= (3, 0):
    str_conversion = lambda x:bytes(x,'utf-8')
else:
    str_conversion = lambda x:x