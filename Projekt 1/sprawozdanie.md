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
