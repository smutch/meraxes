import sys
import os
import socket

if os.path.exists("deps/custom.py"):
    from custom import deps
elif 'coepp' in socket.gethostname():
    from coepp import deps
else:
    if sys.platform == 'darwin':
        from mac import deps
    elif sys.platform == 'linux2':
        from gstar import deps
