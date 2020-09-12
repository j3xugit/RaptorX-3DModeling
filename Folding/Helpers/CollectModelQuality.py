import numpy as np
import sys
import os

"""
deepScoreFile=sys.argv[1]
fh=open(deepScoreFile, 'r')
content=[ line.strip() for line in list(fh) ]
fh.close()
"""

content = list(sys.stdin)

name = content[8].split()[3]
quality = content[13].split()[ -7: ]

line = '\t'.join( [name] + quality )

print line

