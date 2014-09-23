import sys
import os

if os.path.exists("deps/custom.py"):
    from custom import deps
else:
    if sys.platform == 'darwin':
        from mac import deps
    elif sys.platform == 'linux2':
        from gstar import deps
