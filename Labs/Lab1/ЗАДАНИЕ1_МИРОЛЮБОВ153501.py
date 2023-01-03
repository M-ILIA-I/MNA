import numpy

C = numpy.array([[0.2, 0, 0.2, 0, 0],
                [0, 0.2, 0, 0.2, 0],
                [0.2, 0, 0.2, 0, 0.2],
                [0, 0.2, 0, 0.2, 0],
                [0, 0, 0.2, 0, 0.2]])

D = numpy.array([[2.33, 0.81, 0.67, 0.92, -0.53],
                [-0.53, 2.33, 0.81, 0.67, 0.92],
                [0.92, -0.53, 2.33, 0.81, 0.67],
                [0.67, 0.92, -0.53, 2.33, 0.81],
                [0.81, 0.67, 0.92, -0.53, 2.33]])


b = numpy.array([4.2, 4.2, 4.2, 4.2, 4.2])
b = b.transpose()
A = 8 * C + D
X = numpy.zeros(len(A))


#Прямой обход метод Гаусса
def StraightRun(matrix, b):
   for numRow, row in enumerate(matrix):
# где numRow равен номеру строки, а row это сама строка матрицы
       divider = row[numRow]
       row /= divider
       b[numRow] /= divider
       bfactor = b[numRow]
# вычитаем из всех нижележащих строчек
       for lowerRow in range(numRow + 1, len(matrix)):
           factor = matrix[lowerRow][numRow]
# элемент строки в колонке nrow
           matrix[lowerRow] -= factor * row
# вычитаем, чтобы получить ноль в колонке nrow
           b[lowerRow] -= factor * bfactor
   return matrix, b

# Максимальный элемент по столбцу
def straightRunClmn(matrix, b):
    for nrow in range(len(matrix)):  # nrow равен номеру строки

        pivot = nrow + numpy.argmax(abs(matrix[nrow:, nrow]))
        if pivot != nrow:
            matrix[[nrow, pivot]] = matrix[[pivot, nrow]]  # swap
            change = b[pivot]
            b[pivot] = b[nrow]
            b[nrow] = change
        row = matrix[nrow]
        divider = row[nrow]  # диагональный элемент
        if abs(divider) < 1e-10:
            # почти нуль на диагонали.
            # Продолжать не имеет смысла, результат счёта неустойчив
            print(f"Matrix is incompatible. Max element in column {nrow}: {divider:.3g}")
            exit()
        row /= divider  # делим на диагональный элемент.
        b[nrow] /= divider
        bfactor = b[nrow]
        # теперь надо вычесть приведённую строку из всех нижележащих строчек
        for lowerRow in range(nrow + 1, len(matrix)):
            factor = matrix[lowerRow][nrow]
            # элемент строки в колонке nrow
            matrix[lowerRow] -= factor * row
            # вычитаем, чтобы получить ноль в колонке nrow
            b[lowerRow] -= factor * bfactor

#Поиск макс элемента по всей матрице
def straightRunMax(matrix, b):
   for nrow in range(len(matrix)): # nrow равен номеру строки
       mrow, mcol = MaxElement(matrix, nrow)
       if mrow != nrow:
           matrix[[nrow, mrow]] = matrix[[mrow, nrow]] # swap
           change = b[mrow]
           b[mrow] = b[nrow]
           b[nrow] = change
       swapClmn(matrix, nrow, mcol)
       row = matrix[nrow]
       divider = row[nrow]  # диагональный элемент
       if abs(divider) < 1e-10:
           # почти нуль на диагонали. Продолжать не имеет смысла, результат счёта неустойчив
           print(f"Matrix is incompatible. Max element in matrix: {divider:.3g}")
           exit()
       row /= divider # делим на диагональный элемент.
       b[nrow] /= divider
       bfactor = b[nrow]
       # теперь надо вычесть приведённую строку из всех нижележащих строчек
       for lower_row in range(nrow + 1, len(matrix)):
           factor = matrix[lower_row][nrow]
  # элемент строки в колонке nrow
           matrix[lower_row] -= factor * row
# вычитаем, чтобы получить ноль в колонке nrow
           b[lower_row] -= factor * bfactor

def MaxElement(A, k):
   maximum = A[k-1][k-1]
   maxIndex = [k-1, k-1]
   for i in range(k-1, len(A)):
       for j in range(k-1, len(A)):
           if maximum < A[i][j]:
               maximum = A[i][j]
               maxIndex = [i, j]
   return maxIndex

#Перестановка колонок
def swapClmn(a, i, j):
   for k in range(len(a)):
       a[k][i], a[k][j] = a[k][j], a[k][i]

#Проверка на 0 на главной диагонали
def zerOnDiag(matrix, b):
   nstr = len(matrix)
   for i in range(0, nstr):
       if matrix[i][i] == 0:
           check = True
           for j in range(0, nstr):
               if matrix[i][j] != 0 and matrix[j][i] != 0:
                   swapClmn(matrix, i, j)
                   check = False
                   break
           if check:
               if b[i] == 0:
                   print('The system has infinitive amount of solutions')
               else:
                   print('The system has no solutions')
               exit()

# перебор строк в обратном порядке
def revRun(a, b):
   n = len(a)
   for k in reversed(range(0, n)):
       X[k] = (b[k] - sum(a[k][i] * X[i] for i in range(k + 1, n)))/a[k][k]

#Метод Гаусса поиск максимального по столбцу
def GaussMaxcolumn (A, b):
   zerOnDiag(A, b)
   straightRunClmn(A, b)
   revRun(A,b)
   for i in range(len(X)):
       print("%.4f" % X[i], end=" ")
   print('\n')

#Метод Гаусса
def Gauss(A, b):
   zerOnDiag(A, b)
   StraightRun(A, b)
   revRun(A, b)
   for i in range(len(X)):
       print("%.4f" % X[i], end=" ")
   print('\n')

def GaussMax(A, b):
   zerOnDiag(A, b)
   straightRunMax(A, b)
   revRun(A, b)
   for i in range(len(X)):
       print("%.4f" % X[i], end=" ")
   print('\n')


#print(numpy.linalg.solve(A,b))
#print(Gauss(A,b))
print(GaussMax(A,b))
#print(GaussMaxcolumn(A,b))


