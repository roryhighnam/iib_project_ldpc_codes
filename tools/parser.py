import pstats

p = pstats.Stats('out_profile')
p.sort_stats('cumulative').print_stats(10)
