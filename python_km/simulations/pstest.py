import numpy as np
import pointsource

import time

st = time.clock()

ps = pointsource.DiMatteo()
#ps.flux_max = 50.0
psm = ps.generate_map(5.0, 256, 5.0, 256, 120.0, 325.0, 64)

et = time.clock()

print et-st
