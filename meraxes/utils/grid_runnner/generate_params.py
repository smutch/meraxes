import numpy as np

grid_dim = 4
n_param = 5

SfEfficiency = np.logspace(-3, np.log10(0.30), grid_dim)
SnEjectionEff = np.linspace(0.1, 1, grid_dim)
SnReheatEff = np.linspace(1, 20, grid_dim)
ReincorporationEff = np.linspace(0,1,grid_dim)
MergerTimeFactor = np.linspace(0, 3, grid_dim)


with open("grid_params.txt", "w") as fd:
    fd.write("{:d}\n".format(n_param**grid_dim))
    for i in xrange(grid_dim):
        for j in xrange(grid_dim):
            for k in xrange(grid_dim):
                for l in xrange(grid_dim):
                    for m in xrange(grid_dim):
                        fd.write("{:g} {:g} {:g} {:g} {:g}\n"
                                 .format(SfEfficiency[i], SnEjectionEff[j],
                                         SnReheatEff[k], ReincorporationEff[l],
                                         MergerTimeFactor[m]))
