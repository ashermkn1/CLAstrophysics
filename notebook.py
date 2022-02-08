# %%
import matplotlib.pyplot as plt
from utils import *

# %% create dataframe
orbit = vector_ephemeris_to_dataframe("Earth eph.txt")
orbit.columns = [x.lower() for x in orbit.columns]
orbit = generate_velocity(generate_distance(generate_days(orbit)))

# %% distance plot
orbit.plot.scatter(x='days', y='distance')

# %% distance vs velocity plot, shows inverse relationship
orbit.plot.scatter(x='distance', y='velocity')

# %% orbital characteristics
# perihelion/aphelion
peri = orbit['distance'].min()
aph = orbit['distance'].max()
# semi-major axis
a = (peri + aph) / 2
# eccentricity
e = (aph - peri) / (aph + peri)
# semi-minor axis

res = fit_sin(orbit.days, orbit.distance)
period = res['period']

# %% 3d orbital plot
fig = plt.figure(figsize=(16, 9))
ax = plt.axes(projection="3d")
ax.scatter3D(orbit.x, orbit.y, orbit.z, color='blue')
plt.show()
# %% sine fit curve
res = fit_sin(orbit.days, orbit.distance)
fig = plt.figure(figsize=(16, 9))
plt.plot(orbit.days, orbit.distance, label="orbit")
plt.plot(orbit.days, res['fitfunc'](orbit.days), "r--", label="y fit curve")
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

    res = fit_sin(eph.days, eph.distance)
    period = res['period']
    data['period'].append(period)
    data['a'].append(a)

kepler = pd.DataFrame(data)
kepler = kepler.assign(boom=kepler.a.pow(3) / kepler.period.pow(2))
