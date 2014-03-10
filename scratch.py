import math
n = 922
p = 23.
for my_rank in xrange(int(p)):
    print "core ", my_rank, ": start - ", (math.floor(n/p) * my_rank), "end - ", (math.floor((n/p) * (my_rank+1)))
