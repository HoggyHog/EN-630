import numpy as np
import pandas as pd
from tqdm import tqdm

# Given latitude and longitude
latitude = 26.35  # in degrees
longitude = 73.05  # in degrees
standard_time_longitude = 75  # Assume IST (UTC+5:30) -> 75°E

# Arbitrary tilt and surface azimuth angles
beta = 30  # Tilt angle in degrees
gamma = 0  # Surface azimuth angle (0 = facing south)

# Solar constant
Isc = 1367  # W/m^2

# Load hourly data of GHI and DHI (assuming a DataFrame with columns 'day', 'hour', 'GHI', 'DHI')
data = pd.read_csv("630-datafile.csv")  # Replace with actual data file
new_columns=['Hour','GHI','DHI','Temperature','Wind Speed','IT']



new_df= pd.DataFrame(columns=new_columns)

for i in tqdm(range(data.shape[0])):

    entry=data.iloc[i].copy()

    # Compute day of the year
    n = entry['Hour']//24

    # Compute extraterrestrial solar radiation (I'sc)
    Isc_prime = Isc * (1 + 0.033 * np.cos(np.radians(360 * n / 365)))

    # Compute equation of time correction
    B = np.radians((n - 1) * (360 / 365))
    E = 229.18 * (0.000075 + 0.001868 * np.cos(B) - 0.032077 * np.sin(B) -
                0.014615 * np.cos(2 * B) - 0.04089 * np.sin(2 * B))

    # Compute Local Apparent Time (LAT)

    LAT = entry['Hour'] + (4 * (standard_time_longitude - longitude) + E) / 60

    # Compute hour angle (ω)
    omega = (LAT - 12) * 15  # In degrees

    # Compute declination angle (δ)
    delta = np.radians(23.45 * np.sin(np.radians((284 + n) * (360 / 365))))

    # Compute zenith angle (θz)
    phi = np.radians(latitude)
    cos_theta_z = np.cos(delta) * np.cos(np.radians(omega)) * np.cos(phi) + np.sin(delta) * np.sin(phi)
    theta_z = np.degrees(np.arccos(cos_theta_z))

    # Compute beam normal radiation (Ibn)
    Ibn = (entry['GHI'] - entry['DHI']) / np.maximum(cos_theta_z, 1e-6)  # Avoid division by zero

    # Compute beam radiation (Ib)
    Ib = Ibn * cos_theta_z

    # Compute incidence angle (θ)
    cos_theta = (np.sin(np.radians(beta)) * np.cos(np.radians(gamma)) * (np.cos(delta) * np.cos(np.radians(omega)) * np.sin(phi) - np.sin(delta) * np.cos(phi)) +
                np.sin(np.radians(beta)) * np.sin(np.radians(gamma)) * np.cos(delta) * np.sin(np.radians(omega)) +
                np.cos(np.radians(beta)) * (np.cos(delta) * np.cos(np.radians(omega)) * np.cos(phi) + np.sin(delta) * np.sin(phi)))

    # Compute rb factor
    rb = cos_theta / cos_theta_z

    # Compute rd and rr factors
    rd = (1 + np.cos(np.radians(beta))) / 2
    rho = 0.2  # Ground reflectance (typical value)
    rr = rho * (1 - np.cos(np.radians(beta))) / 2

    # Compute IT (Radiation on tilted surface)
    entry['IT'] = Ib * rb + entry['DHI'] * rd + (Ib + entry['DHI']) * rr

    entries_to_append = [entry.to_frame().T] * len(new_columns)

    
    new_df = pd.concat([new_df, entry.to_frame().T], ignore_index=True)

# Save the results to a new CSV file
new_df.to_csv("tilted_radiation.csv", index=False)

print("Calculation complete. Results saved to 'tilted_radiation.csv'")