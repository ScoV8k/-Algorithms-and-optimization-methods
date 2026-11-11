# Projekt 1.

**Autor:** Mikołaj Taudul  
**Nr albumu:** 318735

## Cel projektu

Celem projektu jest wyznaczenie położenia odbiornika GPS na podstawie współrzędnych satelitów i czasów nadejścia ich sygnałów. Model zapisany jest w kartezjańskim układzie współrzędnych i rozwiązany jako zadanie nieliniowej optymalizacji bez ograniczeń metodą najmniejszych kwadratów.

## Dane i założenia

W obliczeniach przyjęto promień Ziemi $R=6\,378\,137\ \mathrm{m}$ i prędkość sygnału $c=299\,792\,458\ \mathrm{m/s}$. Dane wejściowe pochodzą z pliku `dane_35.csv` i zawierają dla każdego satelity szerokość geograficzną $\theta_i$, długość geograficzną $\phi_i$, wysokość $h_i$ oraz czas nadejścia sygnału $t_i$.

Współrzędne kątowe podane w stopniach przeliczamy na radiany. Przyjmujemy geocentryczny układ kartezjański: początek w środku Ziemi, oś $z$ wskazuje biegun północny, osie $x$ i $y$ leżą w płaszczyźnie równika. Pomijamy poprawki zegara satelity i opóźnienia propagacyjne; błąd zegara odbiornika jest niewiadomą modelu.

## 1. Układ równań określający położenie odbiornika

Niech $N$ oznacza liczbę obserwowanych satelitów. Odległość satelity $i$ od środka Ziemi wynosi $r_i=R+h_i$. Jego współrzędne kartezjańskie to:

$$\begin{aligned} X_i &= r_i\cos\theta_i\cos\phi_i,\\ Y_i &= r_i\cos\theta_i\sin\phi_i,\\ Z_i &= r_i\sin\theta_i. \end{aligned}$$

Po podstawieniu danych z `dane_35.csv` otrzymujemy:

| Satelita | $X_i$ [m] | $Y_i$ [m] | $Z_i$ [m] |
|:--|--:|--:|--:|
| G03 | 12414141.480198 | 9135171.535709 | 21437139.669115 |
| G04 | 22643801.681892 | 2451692.264177 | 13848415.405679 |
| G06 | 4116441.926306 | -15252490.747656 | 21391010.798955 |
| G09 | 25276594.087966 | -7084703.769228 | 4244795.492685 |
| G17 | 21376048.112201 | -13426194.374377 | 9087674.073059 |
| G19 | 12161510.938526 | -15830062.166174 | 17202456.657577 |
| G28 | -7092251.398585 | 13623370.945994 | 21667874.056291 |
| G31 | 2518645.208649 | 18986948.032264 | 17970843.387910 |

Szukane położenie odbiornika oznaczamy $(x,y,z)$, a błąd zegara odbiornika $b$. Przyjmując wzór z prezentacji z wykładu, równanie sfery z niedokładnością zegara ma postać:

$$(x - X_i)^2 + (y - Y_i)^2 + (z - Z_i)^2 = [ (\tilde{t}_r + b - t_{nad,i})c ]^2, \qquad i = 1, 2, \dots, N$$

gdzie $\tilde{t}_r$ to czas odbioru, $t_{nad,i}$ to czas nadania przez satelitę $i$, a $b$ to błąd zegara. Czas z danych $t_i = \tilde{t}_r - t_{nad,i}$ jest zmierzonym czasem propagacji, więc po spierwiastkowaniu:

$$\sqrt{(x-X_i)^2+(y-Y_i)^2+(z-Z_i)^2} = c(t_i + b), \qquad i = 1, 2, \ldots, N$$

W pliku mamy $N=8$ obserwacji. Do wyznaczenia czterech niewiadomych potrzeba minimum czterech satelitów, więc układ jest nadokreślony i rozwiążemy go metodą najmniejszych kwadratów.

## 2. Sformułowanie zadania optymalizacji

W rzeczywistych pomiarach równania z pkt. 1 nie muszą być spełnione dokładnie. Definiujemy resztę pomiarową:

$$
f_i(x,y,z,b)=\sqrt{(x-X_i)^2+(y-Y_i)^2+(z-Z_i)^2}-c(t_i+b),
\qquad i=1,\ldots,N.
$$

Wektor zmiennych decyzyjnych to $\mathbf{q}=(x,y,z,b)^T$. Zadanie optymalizacji:

$$
{
\min_{\mathbf{q}\in\mathbb{R}^4}\;\Phi(\mathbf{q})
=\sum_{i=1}^{N} f_i(\mathbf{q})^2
=\sum_{i=1}^{N}\left[\sqrt{(x-X_i)^2+(y-Y_i)^2+(z-Z_i)^2}-c(t_i+b)\right]^2.
}
$$


Mamy $N=8$ reszt i 4 niewiadome. Jest to zadanie optymalizacji nieliniowej bez ograniczeń.

## 3. Wyznaczenie położenia metodą Levenberga-Marquardta

Do rozwiązania użyto Pythona z funkcją `scipy.optimize.least_squares` z opcją `method="lm"` (metoda Levenberga-Marquardta). Dla poprawy skalowania zmienną $b$ zastąpiono zmienną $B=cb$ wyrażoną w metrach, więc reszty mają postać:

$$
f_i(x,y,z,B)=\sqrt{(x-X_i)^2+(y-Y_i)^2+(z-Z_i)^2}-ct_i-B.
$$

Punkt startowy: $(6\,378\,137,\; 0,\; 0,\; 0)$

Solver zwrócił położenie:

$$
(x,y,z)=(3\,662\,131.489,\; 1\,384\,295.286,\; 5\,035\,349.606)\ \mathrm{m}
$$

co odpowiada współrzędnym:

$$
\theta=52.134592°,\qquad \phi=20.706723°,\qquad h=125.0\ \mathrm{m}.
$$

Błąd zegara: $b = 3.69\cdot10^{-14}$ s (praktycznie zero).\
Solver zakończył po 31 ewaluacjach.

## 4. Wpływ parametrów i zaburzeń danych

### 4.1. Punkt startowy

Przetestowano punkty startowe położone w pobliżu Ziemi oraz punkty znacznie od niej oddalone:

| Punkt startowy | $\theta$ [°] | $\phi$ [°] | $h$ [m] | Koszt [m²] | nfev |
|:--|--:|--:|--:|--:|--:|
| Równik $(R,0,0,0)$ | 52.1346 | 20.7067 | 125.0 | $4.08\cdot10^{-9}$ | 31 |
| Biegun N $(0,0,R,0)$ | 52.1346 | 20.7067 | 125.0 | $4.08\cdot10^{-9}$ | 26 |
| Środek Ziemi $(0,0,0,0)$ | 52.1346 | 20.7067 | 125.0 | $4.08\cdot10^{-9}$ | 106 |
| Punkt w przestrzeni $(R,R,R,0)$ | 52.1346 | 20.7067 | 125.0 | $4.08\cdot10^{-9}$ | 26 |
| Skrajny $(10^8, 10^8, 10^8, 0)$ | 59.6919 | 10.3142 | $4.1\cdot10^{11}$ | $1.93\cdot10^{13}$ | 197 |
| Skrajny $(10^9, 10^9, 10^9, 0)$ | 59.6914 | 10.3116 | $2.9\cdot10^{11}$ | $1.93\cdot10^{13}$ | 190 |

Dla punktów startowych leżących w rozsądnej odległości od Ziemi (nawet ze środka Ziemi, gdzie w pierwszym kroku odległości do satelitów są symetryczne, czy z punktu $(R,R,R,0)$) algorytm zbiega do właściwego minimum globalnego z kosztem $4.08\cdot10^{-9}\ \mathrm{m}^2$. 

Przy zadaniu bardzo dużych współrzędnych początkowych ($10^8$ i $10^9\ \mathrm{m}$) metoda Levenberga-Marquardta nie zbiega do poprawnego rozwiązania. Wartość funkcji celu rzędu $1.93\cdot10^{13}\ \mathrm{m}^2$ pokazuje, że krok optymalizacji w tak dużej odległości od sfer pomiarowych staje się nieefektywny i solver zatrzymuje się w obszarze o płaskim gradiencie.

### 4.2. Wpływ kryteriów zatrzymania metody

Zbadano zachowanie solvera przy ograniczeniu maksymalnej liczby wywołań funkcji reszt (`max_nfev`) oraz przy zmianie tolerancji zbieżności. Błąd położenia liczono jako odległość euklidesową od rozwiązania dokładnego (uzyskanego dla tolerancji $10^{-15}$).

**Wpływ limitu wywołań funkcji reszt**

| Limit `max_nfev` | Faktyczne nfev | Koszt [m²] | Błąd położenia [m] | Wysokość [m] |
|--:|--:|--:|--:|--:|
| 10 | 11 | $2.19\cdot10^{11}$ | $4.47\cdot10^5$ | 416166.38 |
| 15 | 16 | $3.29\cdot10^7$ | $4.88\cdot10^3$ | 4642.92 |
| 20 | 21 | $4.67\cdot10^{-1}$ | 0.595 | 125.55 |
| 25 | 26 | $4.08\cdot10^{-9}$ | $1.41\cdot10^{-7}$ | 125.00 |
| 30 | 31 | $4.08\cdot10^{-9}$ | $0.00$ | 125.00 |

**Wpływ tolerancji zbieżności** (`ftol=xtol=gtol`)

| Tolerancja | nfev | Koszt [m²] | Błąd położenia [m] | Wysokość [m] |
|--:|--:|--:|--:|--:|
| $10^{-1}$ | 21 | $4.67\cdot10^{-1}$ | 0.595 | 125.55 |
| $10^{-3}$ | 26 | $4.08\cdot10^{-9}$ | $1.41\cdot10^{-7}$ | 125.00 |
| $10^{-6}$ | 26 | $4.08\cdot10^{-9}$ | $1.41\cdot10^{-7}$ | 125.00 |
| $10^{-10}$ | 31 | $4.08\cdot10^{-9}$ | $0.00$ | 125.00 |
| $10^{-15}$ | 37 | $4.08\cdot10^{-9}$ | $0.00$ | 125.00 |

Przy limicie 20 wywołań błąd wynosi ok. 0,6 m, natomiast 26 ewaluacji wystarcza do osiągnięcia precyzji poniżej mikrometra. Tolerancja na poziomie $10^{-3}$ zapewnia już zbieżność do właściwego rozwiązania; zaostrzanie kryteriów poniżej $10^{-10}$ wymusza jedynie dodatkowe ewaluacje bez zmiany wyznaczonej pozycji.

### 4.3. Wpływ zaburzeń danych

Do czasów nadejścia sygnału dodano losowy szum o rozkładzie normalnym $\mathcal{N}(0, \sigma^2)$. Wyniki pojedynczych realizacji dla $\sigma \in \{1, 10, 100\}\ \mathrm{ns}$ zestawiono w tabeli:

| Odchylenie szumu [ns] | $\theta$ [°] | $\phi$ [°] | $h$ [m] | Koszt [m²] | Błąd położenia [m] | nfev |
|--:|--:|--:|--:|--:|--:|--:|
| 1 | 52.134593 | 20.706720 | 125.61 | $8.38\cdot10^{-2}$ | 0.67 | 31 |
| 10 | 52.134618 | 20.706729 | 121.56 | $1.19\cdot10^{1}$ | 4.51 | 31 |
| 100 | 52.134606 | 20.706922 | 129.84 | $2.79\cdot10^{3}$ | 14.53 | 36 |

Ponieważ pojedyncze losowanie dla małej próby podlega sporym wahaniom losowym, wykonano dodatkowo symulację Monte Carlo (100 powtórzeń) dla $\sigma = 10\ \mathrm{ns}$. Średni błąd położenia wyniósł 5,57 m. Odchylenia standardowe poszczególnych składowych wyniosły: 2,97 m dla szerokości, 2,10 m dla długości oraz 5,19 m dla wysokości. Składowa wysokościowa jest najbardziej wrażliwa na błędy pomiaru czasu, co wynika z geometrii obserwacji.

## 5. Weryfikacja położenia w Google Maps

Wyznaczone współrzędne geograficzne to $\theta = 52°\;8'\;4{,}5''\;\mathrm{N}$, $\phi = 20°\;42'\;24{,}2''\;\mathrm{E}$, a wysokość $h = 125\ \mathrm{m}$ n.p.m.

Lokalizacja na mapie: [Google Maps](https://www.google.com/maps?q=52.134591666,20.706723333)

Wyznaczony punkt znajduje się w miejscowości Brwinów przy ulicy Grodziska 46.


## Wnioski

Metoda Levenberga-Marquardta pozwoliła wyznaczyć położenie odbiornika GPS z normą reszt na poziomie $9\cdot10^{-5}\ \mathrm{m}$ i błędęm zegara wynoszącym $b \approx 0$

Przeprowadzona analiza wrażliwości wykazała, że:
- wybór punktu startowego w otoczeniu Ziemi nie wpływa na końcowy wynik optymalizacji (zmienia się jedynie liczba iteracji), natomiast start z odległości rzędu $10^8\ \mathrm{m}$ powoduje brak zbieżności do właściwego minimum,
- w teście kryteriów stopu 20 wywołań funkcji reszt pozwala określić pozycję z dokładnością do 0,6 m, a 26 wywołań daje błąd poniżej mikrometra; ustawienie tolerancji na poziomie $10^{-3}$ jest w pełni wystarczające,
- błędy pomiaru czasu propagują się na błąd pozycji z wyraźną dominacją błędu wysokościowego (w próbie Monte Carlo dla szumu 10 ns odchylenie standardowe wysokości wyniosło 5,19 m przy ok. 2–3 m dla składowych poziomych).


## Załącznik B. Pełny kod źródłowy

Poniżej zamieszczono pełną treść skryptów wykorzystanych w projekcie. Osobne pliki źródłowe oraz archiwum z kodem i danymi wejściowymi dołączono również do materiałów projektu.

### rozwiazanie_gps.py

```python
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
```

### analiza_wplywu.py

```python
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

```

### rozwiazanie_ampl.py

```python
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
```
