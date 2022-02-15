# %%
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from utils import *

# %% clear plots
plt.rc('text', usetex=True)
plt.cla()
plt.clf()
# %% create dataframe
orbit = vector_ephemeris_to_dataframe("Earth eph.txt")
orbit.columns = [x.lower() for x in orbit.columns]
orbit = generate_velocity(generate_distance(generate_days(orbit)))

# %% distance plot
ax = plt.axes()
ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y / 1e8))
ax.yaxis.set_major_formatter(ticks_y)
plt.plot(orbit.iloc[:500].days, orbit[:500].distance)
ax.set_xlabel('Time (days)')
ax.set_ylabel(r'Distance ($10^8$ km)')
plt.show()
# %% distance vs velocity plot, shows inverse relationship
ax = plt.axes()
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / 1e8))
ax.xaxis.set_major_formatter(ticks_x)
ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y / 1e6))
ax.yaxis.set_major_formatter(ticks_y)
plt.plot(orbit.iloc[:400].distance, orbit[:400].velocity)
ax.set_xlabel(r'Distance ($10^8$ km)')
ax.set_ylabel(r'Velocity ($10^6$ km / day)')
plt.show()

# %% orbital characteristics
# perihelion/aphelion
peri = orbit['distance'].min()
aph = orbit['distance'].max()
# semi-major axis
a = (peri + aph) / 2
# eccentricity
e = (aph - peri) / (aph + peri)
# semi-minor axis
b = a * np.sqrt(1 - np.square(e))

res = fit_sin(orbit.days, orbit.distance)
period = res['period']

# %% 3d orbital plot with ellipse fitting
ax = plt.axes()
ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y / 1e8))
ax.yaxis.set_major_formatter(ticks_y)
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / 1e8))
ax.xaxis.set_major_formatter(ticks_x)
# ticks_z = ticker.FuncFormatter(lambda z, pos: '{0:g}'.format(z / 1e3))
# ax.zaxis.set_major_formatter(ticks_z)
# get sample x and then use ellipse equation for y values
xline = np.linspace(-1.5e8, 1.5e8, 100000)
yline = b * np.sqrt(1 - (np.square(xline) / np.square(a)))
# axis labels
# ax.set_zlabel(r'Z Distance ($10^8$km)')
ax.set_ylabel(r'Y Distance ($10^8$km)')
ax.set_xlabel(r'X Distance ($10^3$km)')
# plot ellipse
ax.plot(xline, yline, color='orange', label=r'$y=b\sqrt{1-\frac{x^2}{a^2}}$')
ax.plot(xline, -yline, color='orange', label=r'$y=-b\sqrt{1-\frac{x^2}{a^2}}$')
# plot orbit
ax.plot(orbit.x, orbit.y, color='blue', label='orbit')
plt.legend(loc="best")
plt.show()
# %% sine fit curve
res = fit_sin(orbit.days, orbit.distance)
ax = plt.axes()
ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y / 1e8))
ax.yaxis.set_major_formatter(ticks_y)
plt.plot(orbit.iloc[:500].days, orbit[:500].distance, label='orbit')
ax.set_xlabel('Time (days)')
ax.set_ylabel(r'Distance ($10^8$ km)')
plt.plot(orbit.iloc[:500].days, res['fitfunc'](orbit.iloc[:500].days), "r--", label="best fit curve")
plt.legend(loc="best")
plt.show()
# %% ratio between period and semi-major axis
orbits = ['Merc eph.txt', 'Venus eph.txt', 'Earth eph.txt', 'Mars eph.txt', 'Jupiter eph.txt', 'Saturn eph.txt',
          'Uranus eph.txt',
          'Neptune eph.txt']
data = {'period': [], 'a': []}
for filename in orbits:
    eph = vector_ephemeris_to_dataframe(filename)
    eph.columns = [x.lower() for x in eph.columns]
    eph = generate_velocity(generate_distance(generate_days(eph)))
    # convert km to au
    eph.distance = eph.distance / 1.496e+8
    peri = eph['distance'].min()
    aph = eph['distance'].max()
    # semi-major axis
    a = (peri + aph) / 2
    # eccentricity
    e = (aph - peri) / (aph + peri)
    print(filename.split()[0], a)
    res = fit_sin(eph.days, eph.distance)
    period = res['period']
    data['period'].append(period)
    data['a'].append(a)
#%%
kepler = pd.DataFrame(data)
kepler = kepler.assign(constant=(kepler.a.pow(3) / kepler.period.pow(2)) * pow(10, 6))
kepler.a = kepler.a.pow(3)
kepler.period = kepler.period.pow(2)
kepler.plot.scatter(x="a", y='period')