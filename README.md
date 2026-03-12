# Método Simplex
### Problema: 
##### En el método simplex, en cada iteración es necesario calcular la solución básica actual y otras cantidades asociadas a la base. Así pues, esto implica resolver repetidamente sistemas de ecuaciones lineales (SEL) de la forma Bx=b, donde B es la matriz base formada por las columnas de A correspondientes a las variables básicas. También, se deben resolver sistemas como Bλ=c para obtener los multiplicadores simplex y calcular los costos reducidos. Si estos sistemas se resolvieran calculando la inversa de B en cada paso, el costo computacional sería alto e innecesario. Por esta razón, se utilizan métodos eficientes para resolver SEL, como la factorización LU junto con sustitución hacia adelante y hacia atrás, que permiten resolver varios sistemas con la misma matriz B de forma rápida después de factorizarla una sola vez. De esta manera, los métodos de SEL son fundamentales para que el simplex pueda ejecutarse de manera eficiente en cada iteración.
### 
### Modelo:
##### Para resolverlo, dividimos la matriz A en una base B (columnas de variables básicas) y una parte no básica. El proceso es cíclico y se apoya en la factorización LUPS, que descompone la base mediante la propiedad PB = LU para evitar el cálculo costoso de la inversa y mejorar la estabilidad numérica.
##### En cada iteración, el primer paso es determinar la solución actual del sistema B * xB = b. Para ello, aprovechamos que L y U son triangulares y aplicamos de forma secuencial las funciones de sustitución:
##### y = forsub(L, b(p))
##### xB = backsub(U, y)
##### Esta propiedad de las matrices triangulares permite que el cálculo de los valores de las variables básicas sea extremadamente rápido y preciso.
##### Posteriormente, evaluamos si la solución es óptima mediante los costos reducidos. Primero obtenemos los multiplicadores del sistema B^T * lambda = c_B. Debido a las propiedades de la transposición, el sistema se resuelve como (U^T * L^T * P) * lambda = c_B, lo que en código se traduce como:
##### w = forsub(U', c_B)
##### v = backsub(L', w)
##### lambda(p) = v
##### Con el vector lambda calculado, derivamos los costos reducidos de todo el sistema mediante la operación:
##### r = c - A^T * lambda
##### Si el valor máximo de r es menor o igual a una tolerancia (tol), el algoritmo termina. Si existe un r_j > 0, esa variable entra a la base y calculamos su dirección de impacto d resolviendo B * d = a_j:
##### y = forsub(L, a_j(p))
##### d = backsub(U, y)
##### Finalmente, el método decide qué variable debe salir mediante el test de la razón mínima:
##### theta = min(xB_i / d_i) para d_i > 0
##### Esta propiedad garantiza que, al aumentar la nueva variable, ninguna de las anteriores se vuelva negativa, manteniendo siempre la factibilidad de la solución hasta alcanzar el óptimo.
###  
### Problemas de prueba:
##### Los problemas de prueba son problemas obtenidos de un curso de Optimización.
##### Se debe llamar a la función siguiendo las siguientes definiciones:
| Variable | Entrada |
|----------|---------|
| A | Matriz de coeficientes de las restricciones |
| c | Vector columna de coeficientes de la función objetivo |
| b | Vector columna del lado derecho de las restricciones |
| Bx | Vector de índices iniciales de las variables básicas |

| Variable | Salida |
|----------|--------|
| s | Vector columna con las variables solución |
| z | Valor óptimo de la función objetivo |
| Bx | Vector de variables que quedaron en la base |