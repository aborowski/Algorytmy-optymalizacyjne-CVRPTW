Kompilacja program�w w standardzie C++11. Przyk�adowo z konsoli:
g++ oszczednosciowyv2.cpp -std=c++11 -o oszczednosciowyv2
g++ tabu.cpp -std=c++11 -o tabu

U�ytkowanie:
Programy wymagaj� plik�w wejsciowych zapisanych w formacie Solomona
http://neo.lcc.uma.es/vrp/vrp-instances/description-for-files-of-solomons-instances/
Przyk�adowe instancje testowe znajduj� si� w paczce, oraz na stronie
http://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-with-time-windows-instances/
Pliki wyj�ciowe tworzone s� w zmodyfikowanym formacie wyj�ciowym Solomona (liczba tras oraz ca�kowita d�ugo�� w pierwszym wierszu)

Uruchamianie program�w:
Programy uruchamiane z linii polece� z 2 argumentami
<nazwa pliku programu> <nazwa pliku wej�ciowego> <opcjonalnie liczba wierzcho�k�w do wczytania>
Po przetworzeniu, plik wyj�ciowy zapisany jest w formacie
res-<liczba wierzcho�k�w lub full>-<nazwa pliku>