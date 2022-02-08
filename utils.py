import numpy as np
import pandas as pd
import scipy.optimize


# courtesy of stack overflow: https://stackoverflow.com/a/42322656
# i understand fourier transforms bc of 3b1b but this is blasphemy at its finest
def fit_sin(time, data):
    time = np.array(time)
    data = np.array(data)

    sampling_frequencies = np.fft.fftfreq(len(time), (time[1] - time[0]))
    Fyy = abs(np.fft.fft(data))
    # excluding the zero frequency "peak", which is related to offset
    guess_freq = abs(sampling_frequencies[np.argmax(Fyy[1:]) + 1])
    guess_amp = np.std(data) * 2 ** 0.5
    guess_offset = np.mean(data)
    guess = np.array([guess_amp, 2 * np.pi * guess_freq, 0, guess_offset])

    def sinfunc(t, A, w, p, c): return A * np.sin(w * t + p) + c

    popt, _ = scipy.optimize.curve_fit(sinfunc, time, data, p0=guess)
    amp, omega, p, c = popt
    f = omega / (2. * np.pi)
    def fitfunc(t): return amp * np.sin(omega * t + p) + c
    return {"amp": amp, "omega": omega, "phase": p, "offset": c, "freq": f, "period": 1. / f, "fitfunc": fitfunc}


def generate_days(df: pd.DataFrame):
    return df.assign(days=(df['date'] - df['date'].iloc[0]).dt.days.astype(np.int64))


def generate_distance(df: pd.DataFrame):
    return df.assign(distance=np.sqrt(df[['x', 'y', 'z']].pow(2).sum(axis=1)))


def generate_velocity(df: pd.DataFrame):
    return df.assign(velocity=np.sqrt(df[['vx', 'vy', 'vz']].pow(2).sum(axis=1)))


def ephemeris(filename):
    data = []
    start = end = 0
    # read lines of file
    with open(filename) as eph:
        contents = eph.readlines()
    for index, line in enumerate(contents):
        # search for start of ephemeris
        if "$$SOE\n" == line:
            # column names
            data.append([title.strip() for title in contents[index - 2].split(",")])
            start = index
        # end of ephemeris
        if "$$EOE\n" == line:
            end = index
    # get all the data
    for data_line in contents[start + 1:end]:
        data.append(data_line.split(","))
    # associate data with column names
    return {
        title: [row[index] for row in data[1:]]
        for index, title in enumerate(data[0])
    }


def vector_ephemeris_to_dataframe(filename):
    df = pd.DataFrame(ephemeris(filename))
    # remove extra newline column and julian days
    del df['']
    del df["JDTDB"]
    # trim AD and BC off of dates
    df = df.rename({"Calendar Date (TDB)": "date"}, axis=1)
    df["date"] = df["date"].str.replace(r'^[^0-9]*', '', regex=True).astype(np.datetime64)

    # convert scientific notation strings to actual numbers (only do it for valid columns)

    # first find which valid columns are in the dataframe
    valid = {'X', 'Y', 'Z', 'VX', 'VZ', 'VY'}.intersection(df.columns)
    # then make all of them numbers
    df[list(valid)] = df[list(valid)].astype(np.float64)
    return df


def observer_table_to_dataframe(filename):
    df = pd.DataFrame(ephemeris(filename))
    # delete empty column
    del df['']
    # simplify names
    df.columns = [x.lower().replace(".", "").removesuffix('_(icrf)') for x in df.columns]
    df = df.rename({"date__(ut)__hr:mn:ss": 'date'}, axis=1)
    # change date column to actual dates
    df['date'] = df['date'].astype(np.datetime64)
    # first find which valid columns are in the dataframe
    valid = {'ra', 'dec'}.intersection(df.columns)
    # then make all of them numbers
    df[list(valid)] = df[list(valid)].astype(np.float64)
    return df
