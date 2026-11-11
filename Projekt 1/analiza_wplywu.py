import csv
import numpy as np
from scipy.optimize import least_squares

c = 299_792_458.0
r_earth = 6_378_137.0

with open("dane_35.csv", encoding="utf-8-sig", newline="") as f:
    reader = csv.DictReader(f, delimiter=";")
    data = list(reader)

lat = np.deg2rad([float(row["Szerokość dziesiętna [°]"]) for row in data])
lon = np.deg2rad([float(row["Długość dziesiętna [°]"]) for row in data])
h = np.array([float(row["Wysokość n.p.m. [m]"]) for row in data])
t = np.array([float(row["Czas nadejścia sygnału GPS [s]"]) for row in data])

r = r_earth + h
sat_x = r * np.cos(lat) * np.cos(lon)
sat_y = r * np.cos(lat) * np.sin(lon)
sat_z = r * np.sin(lat)
satellites = np.column_stack((sat_x, sat_y, sat_z))
measured_dist = c * t

guess_start = [r_earth, 0.0, 0.0, 0.0]

def residuals(params, md_override=None):
    xyz, B = params[:3], params[3]
    md = measured_dist if md_override is None else md_override
    return np.linalg.norm(satellites - xyz, axis=1) - md - B

# Wyznaczenie dokladnego rozwiazania referencyjnego (pelna precyzja float64)
ref_res = least_squares(
    residuals, guess_start, method="lm",
    ftol=1e-15, xtol=1e-15, gtol=1e-15, max_nfev=5000
)
REF_XYZ = ref_res.x[:3]
REF_R = np.linalg.norm(REF_XYZ)
REF_LAT = np.degrees(np.arcsin(REF_XYZ[2] / REF_R))

def solve(guess, ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=5000, md_override=None):
    res = least_squares(
        residuals, guess, args=(md_override,), method="lm",
        ftol=ftol, xtol=xtol, gtol=gtol, max_nfev=max_nfev
    )
    x, y, z = res.x[:3]
    R = np.linalg.norm([x, y, z])
    return {
        "lat": np.degrees(np.arcsin(z / R)),
        "lon": np.degrees(np.arctan2(y, x)),
        "alt": R - r_earth,
        "cost": res.cost,
        "nfev": res.nfev,
        "pos_err": np.linalg.norm(res.x[:3] - REF_XYZ),
    }

# 4.1 Punkt startowy
print("=== 4.1 PUNKT STARTOWY ===\n")

punkty = {
    "(R, 0, 0, 0)": [r_earth, 0, 0, 0],
    "(0, 0, R, 0)": [0, 0, r_earth, 0],
    "(0, 0, 0, 0)": [0, 0, 0, 0],
    "(R, R, R, 0)": [r_earth, r_earth, r_earth, 0],
    "(1e8, 1e8, 1e8, 0)": [1e8, 1e8, 1e8, 0],
    "(1e9, 1e9, 1e9, 0)": [1e9, 1e9, 1e9, 0],
}

for label, g in punkty.items():
    wynik = solve(g)
    print(f"  {label}  lat={wynik['lat']}  lon={wynik['lon']}  "
          f"alt={wynik['alt']}  koszt={wynik['cost']}  nfev={wynik['nfev']}")

# 4.2 Test kryteriow zatrzymania
print("\n=== 4.2 TEST KRYTERIÓW ZATRZYMANIA ===\n")

print("--- 4.2.1 Wpływ limitu wywołań funkcji reszt ---")
for max_n in [10, 15, 20, 25, 30]:
    wynik = solve(guess_start, max_nfev=max_n)
    print(f"  max_nfev={max_n} | nfev={wynik['nfev']} | koszt={wynik['cost']} | "
          f"błąd pozycji={wynik['pos_err']} m | wysokość={wynik['alt']} m")

print("\n--- 4.2.2 Wpływ tolerancji zbieżności ---")
for tol in [1e-1, 1e-3, 1e-6, 1e-10, 1e-15]:
    wynik = solve(guess_start, ftol=tol, xtol=tol, gtol=tol)
    print(f"  tolerancja={tol} | nfev={wynik['nfev']} | koszt={wynik['cost']} | "
          f"błąd pozycji={wynik['pos_err']} m | wysokość={wynik['alt']} m")

# 4.3 Zaburzenia danych
print("\n=== 4.3 ZABURZENIA DANYCH ===\n")

np.random.seed(42)
for sigma_ns in [1.0, 10.0, 100.0]:
    sigma_s = sigma_ns * 1e-9
    t_noisy = t + np.random.normal(0, sigma_s, size=t.shape)
    wynik = solve(guess_start, md_override=c * t_noisy)
    print(f"  sigma={sigma_ns} ns  lat={wynik['lat']}  lon={wynik['lon']}  "
          f"alt={wynik['alt']}  dpos={wynik['pos_err']} m  nfev={wynik['nfev']}")

# Monte Carlo dla sigma = 10 ns
print("\n--- Monte Carlo (100x) sigma=10 ns ---")
mc_errs, mc_lats, mc_lons, mc_alts = [], [], [], []
for _ in range(100):
    t_mc = t + np.random.normal(0, 10e-9, size=t.shape)
    wynik = solve(guess_start, md_override=c * t_mc)
    mc_errs.append(wynik["pos_err"])
    mc_lats.append(wynik["lat"])
    mc_lons.append(wynik["lon"])
    mc_alts.append(wynik["alt"])

mc_errs = np.array(mc_errs)
mc_alts = np.array(mc_alts)
deg_to_m = (np.pi / 180.0) * r_earth

print(f"  sredni blad pozycji: {mc_errs.mean()} m  (std={mc_errs.std()}, max={mc_errs.max()})")
print(f"  std lat: {np.std(mc_lats)} deg  ({np.std(mc_lats) * deg_to_m} m)")
print(f"  std lon: {np.std(mc_lons)} deg  ({np.std(mc_lons) * deg_to_m * np.cos(np.radians(REF_LAT))} m)")
print(f"  std alt: {mc_alts.std()} m")
