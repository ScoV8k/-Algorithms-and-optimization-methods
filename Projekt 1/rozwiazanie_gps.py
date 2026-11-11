import csv
import numpy as np
from scipy.optimize import least_squares

c = 299_792_458.0
r_earth = 6_378_137.0

# wczytanie danych
with open("dane_35.csv", encoding="utf-8-sig", newline="") as f:
    reader = csv.DictReader(f, delimiter=";")
    data = list(reader)

lat = np.deg2rad([float(row["Szerokość dziesiętna [°]"]) for row in data])
lon = np.deg2rad([float(row["Długość dziesiętna [°]"]) for row in data])
h = np.array([float(row["Wysokość n.p.m. [m]"]) for row in data])

# wspolrzedne kartezjanskie satelitow
r = r_earth + h
sat_x = r * np.cos(lat) * np.cos(lon)
sat_y = r * np.cos(lat) * np.sin(lon)
sat_z = r * np.sin(lat)
satellites = np.column_stack((sat_x, sat_y, sat_z))

# pseudoodleglosci
t = np.array([float(row["Czas nadejścia sygnału GPS [s]"]) for row in data])
measured_dist = c * t

def error_func(params):
    xyz = params[:3]
    b_dist = params[3]
    geom_dist = np.linalg.norm(satellites - xyz, axis=1)
    return geom_dist - measured_dist - b_dist

# punkt startowy - na rowniku
guess = [r_earth, 0.0, 0.0, 0.0]

res = least_squares(
    error_func, guess,
    method="lm",
    ftol=1e-12, xtol=1e-12, gtol=1e-12,
    max_nfev=5000,
)

if not res.success:
    raise RuntimeError(f"Optymalizacja nie powiodla sie: {res.message}")

x, y, z = res.x[:3]
b_dist = res.x[3]
b_time = b_dist / c

# przeliczenie na wspolrzedne sferyczne
R = np.linalg.norm([x, y, z])
lat_res = np.degrees(np.arcsin(z / R))
lon_res = np.degrees(np.arctan2(y, x))
alt_res = R - r_earth

print(f"Polozenie (XYZ) [m]: x={x}, y={y}, z={z}")
print(f"Polozenie (sferyczne): lat={lat_res}, lon={lon_res}, alt={alt_res} m")
print(f"Blad zegara [s]: {b_time}")
print(f"Norma reszt [m]: {np.linalg.norm(res.fun)}")
print(f"Koszt (1/2 sumy kw.) [m^2]: {res.cost}")
print(f"Ewaluacje: {res.nfev} | Status: {res.message}")