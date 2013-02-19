import littleworkers

commands = ["./test_read GiggleZ_LR 931 8 16 %03d %03d"%(snap, snap) for snap
            in xrange(117)]

lil = littleworkers.Pool(workers=4)
lil.run(commands)
