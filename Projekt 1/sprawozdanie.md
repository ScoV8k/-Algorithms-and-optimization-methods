# Projekt 1. Wyznaczanie położenia odbiornika GPS

## Cel projektu

Celem projektu jest wyznaczenie położenia odbiornika GPS na podstawie współrzędnych satelitów i czasów nadejścia ich sygnałów. Model zostanie zapisany w kartezjańskim układzie współrzędnych, a następnie rozwiązany jako zadanie nieliniowej optymalizacji bez ograniczeń metodą najmniejszych kwadratów.

## Dane i założenia

W obliczeniach przyjęto promień Ziemi na poziomie morza $R=6\,378\,137\ \mathrm{m}$ oraz prędkość sygnału $c=299\,792\,458\ \mathrm{m/s}$. Dane wejściowe pochodzą z pliku `dane_35.csv` i zawierają dla każdego satelity szerokość geograficzną $\theta_i$, długość geograficzną $\phi_i$, wysokość $h_i$ oraz czas nadejścia sygnału $t_i$.

Współrzędne kątowe podane w stopniach należy przed obliczeniami przeliczyć na radiany. Przyjmujemy geocentryczny układ kartezjański: początek układu znajduje się w środku Ziemi, oś $z$ wskazuje północny biegun, a osie $x$ i $y$ leżą w płaszczyźnie równika. Pomijamy poprawki zegara satelity oraz opóźnienia propagacyjne; błąd zegara odbiornika uwzględniamy jako niewiadomą modelu.

## 1. Układ równań określający położenie odbiornika

Niech $N$ oznacza liczbę obserwowanych satelitów. Dla satelity $i$ jego odległość od środka Ziemi wynosi $r_i=R+h_i$. Na podstawie współrzędnych geograficznych otrzymujemy jego współrzędne kartezjańskie:

$$\begin{aligned} X_i &= r_i\cos\theta_i\cos\phi_i,\\ Y_i &= r_i\cos\theta_i\sin\phi_i,\\ Z_i &= r_i\sin\theta_i. \end{aligned}$$

Po podstawieniu danych z pliku `dane_35.csv` współrzędne kartezjańskie satelitów wynoszą:

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

Wartości podano w metrach i zaokrąglono do $10^{-6}$ m.

Szukane położenie odbiornika oznaczamy przez $(x,y,z)$, a nieznane obciążenie (błąd) zegara odbiornika przez $b$. Geometryczna odległość odbiornika od satelity $i$ wynosi:

$$d_i(x,y,z)=\sqrt{(x-X_i)^2+(y-Y_i)^2+(z-Z_i)^2}.$$

Przyjmując wzór z prezentacji „Optymalizacja nieliniowa bez ograniczeń, cz. 1” (str. 14), równanie sfery uwzględniające niedokładność zegara dla pojedynczego satelity ma postać:

$$(x - X_i)^2 + (y - Y_i)^2 + (z - Z_i)^2 = [ (\tilde{t}_r + b - t_{nad,i})c ]^2, \qquad i = 1, 2, \dots, N$$

gdzie:

- $\tilde{t}_r$ to czas odbioru sygnału przez urządzenie,
- $t_{nad,i}$ to rzeczywisty czas nadania sygnału przez satelitę $i$,
- $b$ to obciążenie (błąd) zegara odbiornika.

Czas podany w danych do zadania (oznaczmy go jako $t_i$) to w rzeczywistości zmierzony przez urządzenie czas propagacji sygnału, czyli różnica między zafałszowanym czasem odbioru a czasem nadania:

$$t_i = \tilde{t}_r - t_{nad,i}.$$

Uwzględniając tę relację, ostateczny model pomiaru (po spierwiastkowaniu obu stron równania) można zapisać w postaci układu z czterema niewiadomymi $(x,y,z,b)$:

$$\sqrt{(x-X_i)^2+(y-Y_i)^2+(z-Z_i)^2} = c(t_i + b), \qquad i = 1, 2, \ldots, N$$

W pliku `dane_35.csv` znajduje się $N=8$ obserwacji. Po przeliczeniu współrzędnych geograficznych i wysokości satelitów na układ kartezjański otrzymujemy poniższy układ równań. Współrzędne satelitów podano w metrach i zaokrąglono do $10^{-6}$ m, natomiast czasy $t_i$ przepisano z pliku danych:

$$
\begin{aligned}
\sqrt{(x-12414141.480198)^2+(y-9135171.535709)^2+(z-21437139.669115)^2} &= 299792458(0.067185843691705+b), && \text{G03},\\
\sqrt{(x-22643801.681892)^2+(y-2451692.264177)^2+(z-13848415.405679)^2} &= 299792458(0.069898455112322+b), && \text{G04},\\
\sqrt{(x-4116441.926306)^2+(y+15252490.747656)^2+(z-21391010.798955)^2} &= 299792458(0.077835357287101+b), && \text{G06},\\
\sqrt{(x-25276594.087966)^2+(y+7084703.769228)^2+(z-4244795.492685)^2} &= 299792458(0.077479832211484+b), && \text{G09},\\
\sqrt{(x-21376048.112201)^2+(y+13426194.374377)^2+(z-9087674.073059)^2} &= 299792458(0.078196045709321+b), && \text{G17},\\
\sqrt{(x-12161510.938526)^2+(y+15830062.166174)^2+(z-17202456.657577)^2} &= 299792458(0.075816121217462+b), && \text{G19},\\
\sqrt{(x+7092251.398585)^2+(y-13623370.945994)^2+(z-21667874.056291)^2} &= 299792458(0.077663334653222+b), && \text{G28},\\
\sqrt{(x-2518645.208649)^2+(y-18986948.032264)^2+(z-17970843.387910)^2} &= 299792458(0.072965035300031+b), && \text{G31}.
\end{aligned}
$$

Do wyznaczenia czterech niewiadomych potrzeba co najmniej czterech satelitów; większa liczba pomiarów prowadzi do układu nadokreślonego, który zostanie rozwiązany w kolejnym punkcie metodą najmniejszych kwadratów.

## 2. Sformułowanie zadania optymalizacji

W rzeczywistych pomiarach równania z punktu 1 nie muszą być spełnione dokładnie. Dla satelity $i$ definiujemy resztę pomiarową jako różnicę między obliczoną odległością geometryczną a odległością wynikającą z pomiaru czasu:

$$
f_i(x,y,z,b)=\sqrt{(x-X_i)^2+(y-Y_i)^2+(z-Z_i)^2}-c(t_i+b),
\qquad i=1,\ldots,N.
$$

Szukane zmienne decyzyjne zapisujemy w wektorze $\mathbf{q}=(x,y,z,b)^T$. Metoda najmniejszych kwadratów polega na minimalizacji sumy kwadratów reszt:

$$
\boxed{
\min_{\mathbf{q}\in\mathbb{R}^4}\;\Phi(\mathbf{q})
=\sum_{i=1}^{N} f_i(\mathbf{q})^2
=\sum_{i=1}^{N}\left[\sqrt{(x-X_i)^2+(y-Y_i)^2+(z-Z_i)^2}-c(t_i+b)\right]^2.
}
$$

Dla danych z pliku `dane_35.csv` mamy $N=8$, zatem funkcja celu jest sumą ośmiu kwadratów reszt, a optymalizowane są cztery niewiadome: trzy współrzędne odbiornika oraz błąd zegara. Jest to zadanie optymalizacji nieliniowej bez ograniczeń: nie nakładamy na zmienne dodatkowych równości ani nierówności.

## 3. Wyznaczenie położenia metodą Levenberga–Marquardta

Do rozwiązania zastosowano język Python oraz funkcję `scipy.optimize.least_squares` z biblioteki SciPy, wywołaną z opcją `method="lm"`, która realizuje metodę Levenberga–Marquardta. Współrzędne odbiornika i odległości wyznaczone z pomiarów są wyrażone w metrach, dlatego dla poprawy skalowania zmienną błędu zegara $b$ zastąpiono równoważną zmienną $B=cb$, również wyrażoną w metrach. Wektor optymalizowany przez solver miał postać $(x,y,z,B)^T$, a reszty obliczano jako

$$
f_i(x,y,z,B)=\sqrt{(x-X_i)^2+(y-Y_i)^2+(z-Z_i)^2}-ct_i-B.
$$

Jako punkt startowy przyjęto $(6\,378\,137,0,0,0)$, a tolerancje `ftol`, `xtol` i `gtol` ustawiono na $10^{-12}$. Maksymalną liczbę ewaluacji funkcji ustawiono na 5000. Nie nakładano ograniczeń na zmienne.

Solver zwrócił następujące położenie odbiornika w geocentrycznym układzie kartezjańskim:

$$
(x,y,z)=(3\,662\,131.488750,\ 1\,384\,295.286493,\ 5\,035\,349.606461)\ \mathrm{m}.
$$

Dla przyjętego modelu sferycznego odpowiada to współrzędnym:

$$
	heta=52.134591666447^\circ,\qquad
\phi=20.706723333264^\circ,\qquad
h=124.999980\ \mathrm{m}.
$$

Wyznaczone obciążenie zegara wynosi $b=B/c=3.68860438427647\cdot10^{-14}\ \mathrm{s}$. Norma wektora reszt wyniosła $9.03588\cdot10^{-5}\ \mathrm{m}$, a wartość funkcji celu w konwencji SciPy (`result.cost`, czyli połowa sumy kwadratów reszt) wyniosła $4.08235718\cdot10^{-9}\ \mathrm{m}^2$. Solver zakończył działanie po 31 ewaluacjach, spełniając warunek `xtol`.

## 4. Wpływ parametrów i zaburzeń danych

### 4.1. Punkt startowy

<!-- Uzupełnić opis badanych punktów startowych i porównanie wyników. -->

### 4.2. Dokładność kryterium stopu

<!-- Uzupełnić opis badanych tolerancji i porównanie wyników. -->

### 4.3. Zaburzenia danych

<!-- Uzupełnić sposób wprowadzenia zaburzeń oraz ich wpływ na wyznaczone położenie. -->

## 5. Weryfikacja położenia w Google Maps

<!-- Uzupełnić: współrzędne geograficzne, link lub zrzut ekranu z mapy oraz krótki komentarz. -->

## Wnioski

<!-- Uzupełnić po wykonaniu obliczeń i analiz. -->

## Załącznik A. Kod

Pełny, uruchamialny kod wykorzystany do obliczeń znajduje się w pliku [rozwiazanie_gps.py](rozwiazanie_gps.py). Skrypt wczytuje dane z `dane_35.csv`, wyznacza współrzędne satelitów, rozwiązuje zadanie metodą Levenberga–Marquardta i wypisuje wynik wraz z diagnostyką solvera. Wymaga bibliotek NumPy i SciPy.
