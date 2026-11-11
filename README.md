# Algorytmy i metody optymalizacji

Repozytorium zawiera projekt wyznaczania położenia odbiornika GPS z wykorzystaniem współrzędnych satelitów i czasów nadejścia ich sygnałów.

## Projekt 1: położenie odbiornika GPS

### Metoda

Położenie odbiornika oraz błąd jego zegara wyznaczane są przez minimalizację sumy kwadratów reszt pomiarowych. Problem nieliniowy rozwiązano metodą Levenberga–Marquardta (`scipy.optimize.least_squares`, `method="lm"`). Dane o ośmiu satelitach znajdują się w pliku `Projekt 1/dane_35.csv`.

### Wynik

Otrzymano położenie około **52.134592° N, 20.706723° E**, na wysokości około **125 m**. Punkt odpowiada okolicom Młochowa w województwie mazowieckim. Dla danych bazowych błąd zegara odbiornika jest praktycznie równy zeru.

Analiza wpływu obejmuje test różnych punktów startowych, dwóch kryteriów zatrzymania solvera oraz zaburzeń czasów pomiarowych. Dla punktu startowego w pobliżu Ziemi solver osiąga rozwiązanie z błędem poniżej mikrometra po około 25–30 ewaluacjach. W symulacji Monte Carlo z szumem czasu o odchyleniu standardowym 10 ns średni błąd położenia wyniósł około 5.57 m.

### Pliki

- `Projekt 1/sprawozdanie.md` — opis modelu, wyników i wniosków.
- `Projekt 1/rozwiazanie_gps.py` — rozwiązanie zadania GPS.
- `Projekt 1/analiza_wplywu.py` — analiza punktu startowego, kryteriów zatrzymania i zaburzeń danych.
- `Projekt 1/dane_35.csv` — dane wejściowe dla satelitów.

### Uruchomienie

Wymagane biblioteki: NumPy i SciPy. Uruchom skrypty z katalogu projektu:

```bash
cd "Projekt 1"
python3 rozwiazanie_gps.py
python3 analiza_wplywu.py
```
