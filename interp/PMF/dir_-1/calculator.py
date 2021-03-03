import geom_calc as gc
import sys

gc.read_file1(sys.argv[1])
x1 = gc.get_distance('C4','BR6')
x2 = gc.get_distance('C4','O9')
print('%8.6f' % (x1 - x2))
