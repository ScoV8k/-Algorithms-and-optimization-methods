import csv
import numpy as np
from amplpy import AMPL, add_to_path
import amplpy.modules

amplpy.modules.load()

c        = 299_792_458.0
r_earth  = 6_378_137.0

with open("dane_35.csv", encoding="utf-8-sig", newline="") as f:
    reader = csv.DictReader(f, delimiter=";")
    data = list(reader)

lat = np.deg2rad([float(row["Szerokość dziesiętna [°]"]) for row in data])
lon = np.deg2rad([float(row["Długość dziesiętna [°]"]) for row in data])
h   = np.array([float(row["Wysokość n.p.m. [m]"])       for row in data])
t   = np.array([float(row["Czas nadejścia sygnału GPS [s]"]) for row in data])

r     = r_earth + h
sat_x = r * np.cos(lat) * np.cos(lon)
sat_y = r * np.cos(lat) * np.sin(lon)
sat_z = r * np.sin(lat)
N     = len(data)

model_code = r"""
param N;
param c;

param x_s {1..N};
param y_s {1..N};
param z_s {1..N};
param t   {1..N};

var x := 6378137.0;
var y := 0.0;
var z := 0.0;
var B := 0.0;

minimize LeastSquares:
    sum {i in 1..N} (
        sqrt((x - x_s[i])^2 + (y - y_s[i])^2 + (z - z_s[i])^2)
        - c * t[i] - B
    )^2;
"""

ampl = AMPL()
ampl.eval(model_code)

ampl.get_parameter("N").set(N)
ampl.get_parameter("c").set(c)

ampl.get_parameter("x_s").set_values({i+1: sat_x[i] for i in range(N)})
ampl.get_parameter("y_s").set_values({i+1: sat_y[i] for i in range(N)})
ampl.get_parameter("z_s").set_values({i+1: sat_z[i] for i in range(N)})
ampl.get_parameter("t").set_values({i+1: t[i]       for i in range(N)})

ampl.set_option("solver", "minos")
ampl.set_option("minos_options",
    "optimality_tolerance=1e-12 "
    "feasibility_tolerance=1e-12 "
    "function_precision=1e-15"
)

ampl.solve()

x_opt = ampl.get_variable("x").value()
y_opt = ampl.get_variable("y").value()
z_opt = ampl.get_variable("z").value()
B_opt = ampl.get_variable("B").value()
cost  = ampl.get_objective("LeastSquares").value()

b_opt = B_opt / c

R       = np.linalg.norm([x_opt, y_opt, z_opt])
lat_res = np.degrees(np.arcsin(z_opt / R))
lon_res = np.degrees(np.arctan2(y_opt, x_opt))
alt_res = R - r_earth

sat = np.column_stack((sat_x, sat_y, sat_z))
xyz = np.array([x_opt, y_opt, z_opt])
residuals = np.linalg.norm(sat - xyz, axis=1) - c * t - B_opt
norm_res  = np.linalg.norm(residuals)

print(f"Położenie (XYZ) [m]: x={x_opt}, y={y_opt}, z={z_opt}")
print(f"Położenie (BLH): lat={lat_res}°, lon={lon_res}°, alt={alt_res} m")
print(f"Błąd zegara [s]: {b_opt}")
print(f"Norma reszt [m]: {norm_res}")
print(f"Funkcja celu [m^2]: {cost}")