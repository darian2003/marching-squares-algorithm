MARCHING SQUARES ALGORITHM

Tema 1 APD - Balagiu Darian - 334CB

Implementarea solutiei multi-threading a algoritmului Marching Sqaures a fost realizata folosind functiile "pthread_create" si "pthread_join" din biblioteca <pthread.h>. Sincronizarea threadurilor a fost realizata prin intermediul barierei si a functiilor aferente din biblioteca "pthread_barrier_mac.c".

Pentru a obtine timpi de executie suficient de buni, threadurile vor paraleliza in thread_function 3 pasi ai algoritmului Marching Squares: 
1. Scalarea imaginii (Optional)
2. Sample the grid
3. March the sqaures

Deoarece threadurile vor folosi si vor modifica aceleeasi date din memorie, am creat structura thread_arg care va fi pasata functiei thread_function la crearea threadului. Structura contine un camp unic, thread_id, si un camp general care include datele partajate de catre toate threadurile. Prin gruparea informatiilor generale intr-o structura separata pe care o pasam ca si referinta evitam alocari de memorie redundante.

Toate alocarile structrilor sunt realizate dinamic in main().
Campul "new_image" va fi setat la NULL daca rescalarea imaginii nu este necesara si ii va fi alocat un spitiu de memorie in caz contrar.

In interiorul rutinei thread_function, zonele de memorie ale imaginii vor fi partitionate prin intermediul variabilelor p, q, start si end care sunt calculate similar ca in laborator. Deoarece fiecare thread va modifica o zona diferita de memorie => nu vor putea exista race conditions in interiorul aceluasi pas din algoritm => singurul mijloc de sincronizare de care avem nevoie este bariera (nu exista zone critice, deci nu avem nevoie de mutex).



