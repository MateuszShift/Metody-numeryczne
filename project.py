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

def zadanieB_gauss():
    start_time = time.time()
    #TODO: implementacja metody Gaussa-Seidela
    end_time = time.time()
    
    

    pass


if __name__ == "__main__":
    A,b,x = zadanieA()
    L,U,D = rozbijMacierz(A)
    zadanieB_jacobi(D,L,U,b,x,1000,A)
    
    
    
    
    
     