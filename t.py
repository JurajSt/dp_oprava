import numpy as np
rand = np.random.RandomState(42)
t = 100 * rand.rand(100)
y = np.sin(2 * np.pi * t) + 0.1 * rand.randn(100)

from astropy.stats import LombScargle
frequency, power = LombScargle(t, y).autopower()

import matplotlib.pyplot as plt
plt.plot(frequency, power)
plt.show()


print t
print y