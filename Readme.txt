Kompilacja programów w standardzie C++11. Przyk³adowo z konsoli:
g++ oszczednosciowyv2.cpp -std=c++11 -o oszczednosciowyv2
g++ tabu.cpp -std=c++11 -o tabu

U¿ytkowanie:
Programy wymagaj¹ plików wejsciowych zapisanych w formacie Solomona
http://neo.lcc.uma.es/vrp/vrp-instances/description-for-files-of-solomons-instances/
Przyk³adowe instancje testowe znajduj¹ siê w paczce, oraz na stronie
http://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-with-time-windows-instances/
Pliki wyjœciowe tworzone s¹ w zmodyfikowanym formacie wyjœciowym Solomona (liczba tras oraz ca³kowita d³ugoœæ w pierwszym wierszu)

Uruchamianie programów:
Programy uruchamiane z linii poleceñ z 2 argumentami
<nazwa pliku programu> <nazwa pliku wejœciowego> <opcjonalnie liczba wierzcho³ków do wczytania>
Po przetworzeniu, plik wyjœciowy zapisany jest w formacie
res-<liczba wierzcho³ków lub full>-<nazwa pliku>