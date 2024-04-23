import math
import time
import matplotlib.pyplot as plt

def stworz_macierz(N,a1,a2,a3):  
    macierz =  A = [[0]*N for _ in range(N)] #tworzenie macierzy NxN wypełnionej zerami
    for i in range(N):
        A[i][i]=a1
        if i >= 1:
            A[i][i-1]=a2
            A[i-1][i]=a2
        if i >= 2:
            A[i][i-2]=a3
            A[i-2][i]=a3
        
    return macierz

def create_vector(N,f):
    b = []
    for i in range(N):
        b.append(math.sin(i*(f+1)))
    return b

def zadanieA():
    N = 955 #indeks 193355 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! do zmiany do 955
    a1 = 8 
    a2 = -1
    a3 = -1
    f = 3

    A = stworz_macierz(N,a1,a2,a3)
    b = create_vector(N,f)
    x = [1]*N
    return A,b,x
    
def rozbijMacierz(A):
   
    L = [[0]*len(A) for _ in range(len(A))]
    for i in range(len(A)):
       for j in range(i):
           L[i][j] = A[i][j]
    
    U = [[0]*len(A) for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(i+1,len(A)):
            U[i][j] = A[i][j]
    
    D = [[0]*len(A) for _ in range(len(A))]
    for i in range(len(A)):
        D[i][i] = A[i][i]
    return L,U,D 
    
def dodaj_macierze(macierz1, macierz2):
    macierz = [[0]*len(macierz1) for _ in range(len(macierz1))] #ala zeros()
    for i in range(len(macierz1)):
        for j in range(len(macierz1)):
            macierz[i][j] = macierz1[i][j] + macierz2[i][j]
    return macierz

def odwroc_macierz_diagonalna(macierz):
    
    macierz_inv = [[0]*len(macierz) for _ in range(len(macierz))]
    for i in range(len(macierz)):
        macierz_inv[i][i] = 1/macierz[i][i]
    return macierz_inv

def mnozenie_macierz_wektor(macierz,wektor):
    
    wynikowy_wektor = [0]*len(wektor)
    for i in range(len(macierz)):
        for j in range(len(wektor)):
            wynikowy_wektor[i] += macierz[i][j]*wektor[j]
    return wynikowy_wektor
    
def mnozenie_macierz_wartosc(macierz,wartosc):
        
        wynikowa_macierz = [[0]*len(macierz) for _ in range(len(macierz))]
        for i in range(len(macierz)):
            for j in range(len(macierz)):
                wynikowa_macierz[i][j] = macierz[i][j]*wartosc
        return wynikowa_macierz

def mnozenie_macierzy(macierz1,macierz2):

    wynikowa_macierz = [[0]*len(macierz1) for _ in range(len(macierz1))]
    for i in range(len(macierz1)):
        print(i)
        for j in range(len(macierz2[0])):
            for k in range(len(macierz1)):
                wynikowa_macierz[i][j] += macierz1[i][k]*macierz2[k][j]
    return wynikowa_macierz

def norma_wektora(wektor):
    suma = 0
    for i in range(len(wektor)):
        suma += wektor[i]**2
    return math.sqrt(suma)

def dodaj_wektory(wektor1,wektor2):
    wynikowy_wektor = [0]*len(wektor1)
    for i in range(len(wektor1)):
        wynikowy_wektor[i] = wektor1[i] + wektor2[i]
    return wynikowy_wektor

def odejmij_wektory(wektor1,wektor2):
    wynikowy_wektor = [0]*len(wektor1)
    for i in range(len(wektor1)):
        wynikowy_wektor[i] = wektor1[i] - wektor2[i]
    return wynikowy_wektor

def mnozenie_wektor_wartosc(wektor,wartosc):
    wynikowy_wektor = [0]*len(wektor)
    for i in range(len(wektor)):
        wynikowy_wektor[i] = wektor[i]*wartosc
    return wynikowy_wektor

def metoda_jacobiego(D,L,U,b,x, iteracje, A):
    tab_blad = []
    iter = []
    for i in range(iteracje):
        L_U = dodaj_macierze(L,U) # macierz L+U
        d_inv = odwroc_macierz_diagonalna(D)
        tt = mnozenie_macierz_wektor(L_U,x)
        d_inv_minus = mnozenie_macierz_wartosc(d_inv,-1)
        Mr = mnozenie_macierz_wektor(d_inv_minus,tt) #to jest -D^-1*(L+U)*x(k)
        b_mj = mnozenie_macierz_wektor(d_inv,b) #to jest D^-1*b
        x_nowy = dodaj_wektory(Mr,b_mj) #to jest x(k+1)
        blad_Ar = mnozenie_macierz_wektor(A,x_nowy)
        temp = odejmij_wektory(blad_Ar,b)
        blad = norma_wektora(temp)        
        
        if blad < 1e-9:
            x = x_nowy
            break
        x = x_nowy 
        iter.append(i)
        tab_blad.append(blad)
    
    return i, blad, tab_blad, iter
        
def zadanieB_jacobi(D,L,U,b,x,iter,A):
    start_time = time.time()
    #print("rrr")
    iteracje, blad,tab_blad,iter = metoda_jacobiego(D,L,U,b,x, 1000,A)
    end_time = time.time()
    
    print("Iteracje: ", iteracje)
    print("Blad: ", blad)
    print("Czas: ", end_time - start_time)

    plt.plot(iter,tab_blad)
    plt.title('Metoda Jacobiego')
    plt.yscale('log')
    plt.xlabel('iteracje metody Jacobiego')
    plt.ylabel('blad metody Jacobiego')
    plt.show()

def podstawienie_w_przod(L,b):
    n = len(L)
    y = [0]*n
    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= L[i][j]*y[j]
        y[i] = y[i]/L[i][i]
    return y

def mnozenie_wektor_wartosc(wektor,wartosc):
    wynikowy_wektor = [0]*len(wektor)
    for i in range(len(wektor)):
        wynikowy_wektor[i] = wektor[i]*wartosc
    return wynikowy_wektor

def metoda_gaussa_seidela(D,L,U,b,x,iteracje,A):
    tab_blad = []
    iter = []
    for i in range(iteracje):
        r = mnozenie_macierz_wektor(U,x)
        M = podstawienie_w_przod(dodaj_macierze(D,L),r) 
        M_gs = mnozenie_wektor_wartosc(M,-1)
        b_mgs = podstawienie_w_przod(dodaj_macierze(D,L),b)
        
        x_nowy = dodaj_wektory(M_gs,b_mgs)
        
        temp = odejmij_wektory(mnozenie_macierz_wektor(A,x_nowy),b)
        blad = norma_wektora(temp)
        
        if blad < 1e-9:
            x = x_nowy
            break
        x = x_nowy
        iter.append(i)
        tab_blad.append(blad)
        
    return i, blad, tab_blad, iter

def zadanieB_gauss(D,L,U,b,x,iteracje,A):
    iter = []
    tab_blad = []
    start_time = time.time()
    #TODO: implementacja metody Gaussa-Seidela
    
    i, blad, tab_blad, iter = metoda_gaussa_seidela(D,L,U,b,x,iteracje,A)
        
    end_time = time.time()
    
    print ("Iteracje: ", i)
    print ("Blad: ", blad)
    print ("Czas: ", end_time - start_time)
    
    plt.plot(iter,tab_blad)
    plt.title('Metoda Gaussa-Seidela')
    plt.yscale('log')
    plt.xlabel('iteracje metody Gaussa-Seidela')
    plt.ylabel('blad metody Gaussa-Seidela')
    plt.show()
    
    pass

def podstawienie_wstecz(U,y):
    n = len(U)
    x = [0]*n
    for i in range(n-1,-1,-1): # od n-1 do 0
        suma = 0
        for j in range(i+1,n):
            suma += U[i][j]*x[j]
        x[i] = (y[i] - suma)/U[i][i]
    return x
        
def zadanieC():
    N = 955
    a1 = 3
    a2 = -1
    a3 = -1
    f = 3
    A = stworz_macierz(N,a1,a2,a3)
    b = create_vector(N,f)
    L,U,D = rozbijMacierz(A)
    x = [1]*N
    iter_jacob, blad_jacobi, tab_blad_jacobi, iter_jacobi = metoda_jacobiego(D,L,U,b,x, 20,A)
    iteracje_gs, blad_gs, tab_blad_gs, iter_gs = metoda_gaussa_seidela(D,L,U,b,x,20,A)
    
    plt.plot(iter_jacobi, tab_blad_jacobi, label='Metoda Jacobiego')
    plt.plot(iter_gs, tab_blad_gs, label='Metoda Gaussa-Seidela')
    plt.title('Zmiana błędu w kolejnych iteracjach')
    plt.yscale('log')
    plt.xlabel('Iteracje')
    plt.ylabel('Błąd')
    plt.legend()
    plt.show()

def faktoryzacja_lu(N):
    
    a1 = 3
    a2 = -1
    a3 = -1
    f = 3
    A = stworz_macierz(N,a1,a2,a3)
    b = create_vector(N,f)
    
    U = [row[:] for row in A]
    L = [[0]*N for _ in range(N)]
    
    for i in range(N):
        L[i][i] = 1
    
    for i in range(1, N):
        for j in range(i):
            L[i][j] = U[i][j] / U[j][j]
            for k in range(j, N):  
                U[i][k] -= L[i][j] * U[j][k]

    return L,U,A,b
    
def zadanieD():
    
    start_time = time.time()
    N = 955
    L,U,A,b = faktoryzacja_lu(N)
    
    y = podstawienie_w_przod(L,b)
    x = podstawienie_wstecz(U,y)
    
    residiuum = norma_wektora(odejmij_wektory(mnozenie_macierz_wektor(A,x),b))
    
    end_time = time.time()
    print("Residiuum: ", residiuum)
    print("Czas: ", end_time - start_time)
    
def metoda_faktoryzacji_cala(N):
    L,U,A,b = faktoryzacja_lu(N)
    
    y = podstawienie_w_przod(L,b)
    x = podstawienie_wstecz(U,y)
    
    residiuum = norma_wektora(odejmij_wektory(mnozenie_macierz_wektor(A,x),b))
    
#TODO: wykres czasu od znalezenia rozwiązania w zaleznosci od metody dla macierzy N: 100,500,1000,2000,3000 dla macierzy z zadania A
def stworz_dane_wejsciowe(N):
    a1 = 8
    a2 = -1
    a3 = -1
    f = 3
    A = stworz_macierz(N,a1,a2,a3)
    b = create_vector(N,f)
    x = [1]*N
    return A,b,x

def zadanieE():
    czasy_jacobi = []
    czasy_gs = []
    czasy_faktoryzacja = []
    
    tablica_itr = [100, 500, 1000, 2000, 3000]
    
    for i in range(len(tablica_itr)):
        A,b,x = stworz_dane_wejsciowe(tablica_itr[i])
        L,U,D = rozbijMacierz(A)
        time_jacobi = time.time()
        metoda_jacobiego(D,L,U,b,x, 1000,A)
        time_jacobi_finish = time.time()
        czasy_jacobi.append(time_jacobi_finish - time_jacobi)
        
        time_gs = time.time()
        metoda_gaussa_seidela(D,L,U,b,x,1000,A)
        time_gs_finish = time.time()
        czasy_gs.append(time_gs_finish - time_gs)
        
        time_faktoryzacja = time.time()
        metoda_faktoryzacji_cala(tablica_itr[i])
        time_faktoryzacja_finish = time.time()
        czasy_faktoryzacja.append(time_faktoryzacja_finish - time_faktoryzacja)

    
    return czasy_jacobi, czasy_gs, czasy_faktoryzacja

if __name__ == "__main__":
    A,b,x = zadanieA()
    L,U,D = rozbijMacierz(A)
    #zadanieB_jacobi(D,L,U,b,x,1000,A)
    #zadanieB_gauss(D,L,U,b,x,1000,A)
    #zadanieC()
    #zadanieD()
    czas_jacobi, czas_gs, czas_facto =  zadanieE()
    
    print(czas_jacobi)
    print(czas_gs)
    print(czas_facto)
    
    plt.plot([100,500,1000,2000,3000],czas_jacobi, label='Metoda Jacobiego')
    plt.plot([100,500,1000,2000,3000],czas_gs, label='Metoda Gaussa-Seidela')
    plt.plot([100,500,1000,2000,3000],czas_facto, label='Metoda faktoryzacji')
    plt.title('Czas wykonania metody w zależności od N')
    plt.xlabel('N')
    plt.ylabel('Czas')
    plt.legend()
    plt.show()
    
    
    
    
     