import pstats
import sys
# from pstats import SortKey
p = pstats.Stats(sys.argv[1])
#p.strip_dirs().sort_stats(1).print_stats()
p.sort_stats(1).print_stats()
