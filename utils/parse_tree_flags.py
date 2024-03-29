"""Get the on switches for a merger tree flag integer.

Usage: parse_tree_flags.py <flag>
"""

from docopt import docopt

__author__ = "Simon Mutch"
__date__ = "21 Jan 2014"


class tree_flags(object):

    def __init__(self):
        self.flags = {}
        with open("../src/tree_flags.h", "r") as fd:
            for line in fd:
                line = line.split()
                if len(line) == 0:
                    continue
                if not line[0].startswith('//'):
                    if line[2] in self.flags.keys():
                        self.flags[line[1]] = self.flags[line[2]]
                    else:
                        self.flags[line[1]] = int(line[2])
        self.flags = dict((k, v) for k, v in self.flags.items() if 'TTTP' not in k)

    def parse(self, num):
        print(("Flag {:d} matches:".format(num)))
        for s, v in self.flags.items():
            if (num & v) == v:
                print(s)


if __name__ == '__main__':
    args = docopt(__doc__)
    flags = tree_flags()
    flags.parse(int(args["<flag>"]))
