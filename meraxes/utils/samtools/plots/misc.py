import os
observations_dir = os.path.join(os.path.dirname(os.path.realpath( __file__ )), "observations")

def set_custom_rcParams(rcp):
    rcp['lines.linewidth'] = 2.0
    rcp['axes.color_cycle'] = ['348ABD', '7A68A6', 'A60628', '467821', 'CF4457', '188487', 'E24A33']
