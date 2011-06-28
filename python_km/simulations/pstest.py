import numpy as np
import pointsource

import time

st = time.clock()

ps = pointsource.DiMatteo()
ps.x_num, ps.y_num = (256, 256)
ps.nu_lower, ps.nu_upper, ps.nu_num = (120.0, 325.0, 64)
psm = ps.getfield()

et = time.clock()

print et-st
