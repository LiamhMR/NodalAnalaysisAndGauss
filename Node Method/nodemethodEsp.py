# -*- coding: utf-8 -*-
"""NodalMethod.ipynb

#Eliminación Gaussiana para realizar análisis nodal de un circuito
Método de nodos para análisis de circuitos y uso de eliminación Gaussiana para resolver el sistema de ecuaciones.

El código no ha sido optimizado pero puede usarse de referencia, el alcance de este análisis por nodos ha sido probado 
con los circuitos de ejemplo y puede ser funcional para casos con resistencia, fuentes de voltaje y fuentes
de corriente(siempre que sean independientes) y super nodos. Sin embargo quizá exista discrepancia para circuitos más complejos; 
por ejemplo casos donde haya supernodos que incluyan más de 2 nodos.
Úsese como referencia y no como un solucionador de análisis de nodo.


By
Leandro Maureira.

# Análisis y obtención de nodos
"""

import numpy as np
import json

#Definición de las clases a utilizar
class interlmnt:
    def __init__(self,neighbor,lmntIndex):
      self.name=neighbor #nombre del nodo vecino
      self.lmntIndex=lmntIndex #Indice del elemento intermedio

class nodo:
    def __init__(self,name,con,ref):
      self.name=name
      self.con=con #Arreglo con todos los elementos a los que está enlazados por el nodo
      self.ref=ref #Boolean, true si el nodo es nodo de referencia
      self.way=[] #Arreglo que contiene si el sentido de la corriente entra o sale del nodo en ese sentido. #+1 entra #-1
      self.volt=None #Voltaje del nodo, None si no está definida/resuelta
      self.solve=False #Indica si el nodo está resuelto
      self.neighbor=[] #Lista de vecinos en tuplas de clase interlmnt
      self.supernode=None #Si el nodo es parte de un super nodo indica el nombre del nodo con el que es supernodo
      self.check=False #Indica si el nodo ha sido revisado (Para el caso de supernodos), se inicializa en False por defecto

class lmnt:
    def __init__(self,ind,tpe,val,way):
      self.ind=ind #Indice o numero del elemento

      #Valores para tpe
      #0 fuente de tension dependiente #No implementado
      #1 fuente de tension independiente #Funcional
      #2 resistencias #Funcional
      #3 fuentes dependientes de corriente #No implementado
      #4 fuentes independientes de corriente #Funcional
      self.tpe=tpe

      self.val=val #valor

      self.way=way #Nodo al que apunta la corriente/Nodo al que apunta el negativo

#JSON de ejemplo para modelar un circuito
exampleA= """
{
    "lmnts":[
        [0,1,12,"n1"],
        [1,2,4,0],
        [2,2,6,0],
        [3,2,3,0]
    ],

    "nodes":[
        {
            "name":"n1",
            "con":[0,1],
            "ref":false
        },
        {
            "name":"n2",
            "con":[0,2,3],
            "ref":true
        },
        {
            "name":"n3",
            "con":[1,2,3],
            "ref":false
        }]
}
"""
exampleB="""
{
    "lmnts":[
        [0,1,10,"D"],
        [1,2,3,0],
        [2,2,2,0],
        [3,2,2,0],
        [4,2,5,0],
        [5,2,3,0],
        [6,2,2,0]
    ],

    "nodes":[
        {
            "name":"A",
            "con":[0,1,2],
            "ref":false
        },
        {
            "name":"B",
            "con":[1,3,4],
            "ref":false
        },
        {
            "name":"C",
            "con":[3,5,6],
            "ref":false
        },
        {
            "name":"D",
            "con":[0,4,5],
            "ref":true
        },
        {
            "name":"E",
            "con":[2,6],
            "ref":false
        }
    ]
}
"""

exampleC="""
{
    "lmnts":[
        [0,4,5,"C"],
        [1,2,4,0],
        [2,2,4,0],
        [3,4,2,"B"],
        [4,2,3,0],
        [5,4,4,"B"]
    ],

    "nodes":[
        {
            "name":"A",
            "con":[1,2,3],
            "ref":false
        },
        {
            "name":"B",
            "con":[3,4,5],
            "ref":false
        },
        {
            "name":"C",
            "con":[0,2,4,5],
            "ref":true
        },
        {
            "name":"D",
            "con":[0,1],
            "ref":false
        }
    ]
}
"""

exampleD="""
{
    "lmnts":[
        [0,4,2,"A"],
        [1,2,2,0],
        [2,2,4,0],
        [3,4,7,"C"],
        [4,1,2,"A"],
        [5,2,10,0]
    ],

    "nodes":[
        {
            "name":"A",
            "con":[0,1,4,5],
            "ref":false
        },
        {
            "name":"B",
            "con":[2,3,4,5],
            "ref":false
        },
        {
            "name":"C",
            "con":[0,1,2,3],
            "ref":true
        }]
}
"""



datos=json.loads(exampleD)

#Leer el JSON y almacenar sus datos en listas
#retorna dos listas, una con elementos(fuentes, resistencias,etc) y otra con los puntos de conexión
def modCirc(jsonfile):
  listaElementos=[]
  listaNodos=[]
  for i in datos["lmnts"]:
    elemento = lmnt(i[0],i[1],i[2],i[3])
    listaElementos.append(elemento)
    #print(elemento.val)
  for j in datos["nodes"]:
    nod = nodo(j["name"],j["con"],j["ref"])
    listaNodos.append(nod)

  return listaElementos,listaNodos

elementos,nodos=modCirc(datos)
#Actualiza el valor de la multiplicación que acompaña la variable
def updateVars(listVarName,listVarVal,VarName,VarVal):
  index=0
  for name in listVarName:
    if name==VarName:
      listVarVal[index]+=VarVal
    index+=1
  return listVarVal

#Revisa si un valor se encuentra en listVarNames
def varDefPre(listVarNames,name):
  for n in listVarNames:
    if n==name:
      return True
      break
  return False

#Comprueba si el nodo de nombre name en listaNodos corresponde a un nodo de referencia
def isRefNode(name,listaNodos):
  for n in listaNodos:
    if n.name==name and n.ref==True:
      return True
  return False

#Arma la parteA de la matriz en el orden de ingreso de la ecuación; recibe la lista de nodos, la lista de variables, la lista de valores, y la parteA actual
def updateParteA(nodos,listVarNames,listVarVals,parteA):
  for N in nodos:
        indiceVnn=0
        for Vn in listVarNames:
          #Si nuestro nodo es el mismo que el de la variable
          if N.name==Vn:
            parteA.append(listVarVals[indiceVnn])
          indiceVnn+=1
  return parteA

#Buscada el nodo looking en nodos y devuelve el nodo en cuestion
def getNode(nodos,looking):
  for n in nodos:
    if looking==n.name:
      return n
  return None

#Obtiene las ecuaciones de voltaje del circuito y las devuelve en dos matrices a,b.
#Retorna matriz A,matriz B,lista de nodos actualizada con los valores de nodos resueltos(referencias y voltajes a tierra)
 #y la lista de variables en orden.
#Función válida para fuentes de voltaje, fuentes de corriente y resistencia, no probada en supernodos de más de dos nodos
def getEcx(elementos,nodos):
  #Analisis por nodo, recorriendo cada uno

  #Buscar nodos vecinos
  for CurrentN in nodos:
    #Se resuelve el nodo de referencia
    if CurrentN.ref==True:
      CurrentN.volt=0
      CurrentN.solve=True
    for N in nodos:
      #Compruebo que no comparo al nodo con si mismo
      if CurrentN.name!=N.name:
        #Recorro las conexiones de cada nodo
        for currentCon in CurrentN.con:
          for Con in N.con:
            #Si tienen una conexión en común es un nodo vecino
            if currentCon==Con:
              #Guardamos el nodo vecino junto con su elemento intermedio
              CurrentN.neighbor.append(interlmnt(N.name,currentCon))

  listVarNames=[] #Lista con los nombres de las variables a calcular
  listVarVals=[] #Lista con los valores matriciales de cada variable

  for nod in nodos:
    #Para los nodos que no son el de referencia
    if nod.ref!=True:
      #Se comprueba si existen fuentes de tensión presentes en el nodo
      #Recorremos los vecinos del nodo
      for nb in nod.neighbor:
        #Si existen fuentes de tension el voltaje del nodo es igual al valor del voltaje de la tensión siempre y cuando el voltaje esté a tierra.
        if elementos[nb.lmntIndex].tpe==1 and isRefNode(nb.name,nodos)==True:
          #Sentido de la fuente
          if elementos[nb.lmntIndex].way==nb.name:
            nod.volt=elementos[nb.lmntIndex].val
          else:
            nod.volt=-elementos[nb.lmntIndex].val
          #Marco que esta variable ya está resuelta
          nod.solve=True
          #Si se resolvió el nodo después de haber sido agregado este a la lista de variables la quitaremos de la lista
          listVarNames.remove(nod.name)
          #Si el nodo está resuelto continuo con el siguiente
          break
        #Si tengo una fuente de voltaje no conectada a tierra se es supernodo con el nodo con el que comparte la fuente.(Comprobamos que no hayamos resuelto el nodo antes)
        elif elementos[nb.lmntIndex].tpe==1 and isRefNode(nb.name,nodos)==False and nod.solve!=True:
          nod.supernode=nb.name
        #Si la variable no tiene solución directa el nodo queda como incógnita y pasa a ser una variable de la ecuación
        #Comprobamos que la variable no haya sido agregada previamente antes de agregarla ni que haya sido resuelta ya
        if varDefPre(listVarNames,nod.name)==False and nod.solve!=True:
          listVarNames.append(nod.name)
          listVarVals.append(0)
  print(listVarNames)
  #Reiniciamos el bucle con aquellas variables que ya tengamos resueltas y armamos las matrices
  matrixA=[]
  matrixB=[]
  #Análisis nodo a nodo sin contar los nodos resueltos ni los supernodos
  for nod in nodos:
    #Reiniciamos nuestros valores de las incognitas para la ecuación de nuestro nodo actual
    for index in range(len(listVarVals)):
      listVarVals[index]=0

    if nod.solve!=True and nod.supernode==None:
      parteA=[]
      #Si no existen tensiones se calculan las salidas de corriente=0
      #Armado de ecuacion usando los nodos vecinos
      #Divisiones de resistencias
      nVcn=nod.name #Nombre variable del nodo actual
      Vcn=0 #Variable del nodo actual
      #Divisiones de nodos vecinos
      intVal=0
      for lmntIndex in nod.con:
        if elementos[lmntIndex].tpe==2:
          Vcn+=1/elementos[lmntIndex].val

      #Buscar el nodo vecino; exceptuando super nodos
      if nod.supernode==None:
        for nbor in nod.neighbor:
          Vnn=nod.name
          Vvn=0
          #Para cada nodo
          #Recorremos todos nuestros nodos
          for findN in nodos:
            #Para cada nodo definimos una lista de elementos utilizados
            #Y buscamos el nodo que corresponde a nuestro vecino
            if nbor.name==findN.name:
              lmntIndex=nbor.lmntIndex
              #Si el nodo vecino no está resuelto es una variable de la ecuación
              if findN.solve!=True:
                #Si el elemento entre nuestro nodo actual y el nodo vecino es una resistencia
                if elementos[lmntIndex].tpe==2:
                  Vnn=nbor.name #Nombre de la variable corresponderá al nodo vecino
                  Vvn=-1/elementos[lmntIndex].val #Valor por el que se multiplica nuestro vecino 1/El elemento entre nuestro nodo y el vecino
                  #Actualizamos el valor en la ecuación
                  listVarVals=updateVars(listVarNames,listVarVals,Vnn,Vvn)
              #Si el nodo vecino está resuelto pero el elemento entre estos es una fuente de corriente
              else:
                #Si el nodo vecino está resuelto, pero el elemento entre este y el nodo actual es una fuente de voltaje
                intVal+=-findN.volt/elementos[lmntIndex].val
              #Tanto si el nodo vecino está resuelto o no, si tiene un elemento que es fuente de corriente entre estos entonces:
              if elementos[lmntIndex].tpe==4:
                #Si la corriente indepente apunta anuestro nodo actual la corriente se resta, del caso contrario esta se suma.
                if elementos[lmntIndex].way==nod.name:
                  intVal+=-elementos[lmntIndex].val
                else:
                  intVal+=elementos[lmntIndex].val

      #Actualizamos el valor de la ecuación
      listVarVals=updateVars(listVarNames,listVarVals,nVcn,Vcn)

      #La ecuación final corresponde a A*Vcn+ [for n in Vnn.lenght(Vnn[n]*Vvn[n])] = -intVal
      #La forma de armar la matriz será en el mismo orden en que se anotaron los elementos de entrada pero omitiendo las variables resueltas
      #For para respetar el orden de ingreso y tener una matriz final coherente
      for N in nodos:
        indiceVnn=0
        for Vn in listVarNames:
          #Si nuestro nodo es el mismo que el de la variable
          if N.name==Vn:
            parteA.append(listVarVals[indiceVnn])
          indiceVnn+=1
      #Al final del proceso agregamos la ecuación en forma matricial por partes
      matrixA.append(parteA)
      #Agregamos la igualdad
      matrixB.append(-intVal)

  #Casos para supernodos
  #Primera parte
  #Para este caso se resultan dos ecuaciones, una es que la resta de un nodo con el otro debe ser igual al voltaje.
  for nod in nodos:
    #Si el nodo es parte de un super nodo y no ha sido revisado previamente
    if nod.supernode!=None and nod.check==False:
      #Reiniciamos nuestros valores de las incognitas para la ecuación de nuestro nodo actual
      for index in range(len(listVarVals)):
        listVarVals[index]=0
      parteA=[]
      #Buscamos el vecino que corresponde al supernodo
      for nbor in nod.neighbor:
        #Explicitamos que estemos tratando con el elementos correcto(fuente de voltaje)
        if nod.supernode==nbor.name and elementos[nbor.lmntIndex].tpe==1:
          #Este es nuestro super nodo, nod + nod.supernode

          #Definimos la primera ecuación, nod-nbor=Voltaje del elemento intermedio, es decir (1  -1) (Vx), con el negativo en el nodo al que no apunta el negativo de la fuente
          #Si el negativo de la fuente apunta al nodo actual este es positivo en la matriz
          if elementos[nbor.lmntIndex].way==nod.name:
            listVarVals=updateVars(listVarNames,listVarVals,nod.name,1)
            listVarVals=updateVars(listVarNames,listVarVals,nbor.name,-1)
          #Sino el negativo es el opuesto
          else:
            listVarVals=updateVars(listVarNames,listVarVals,nod.name,-1)
            listVarVals=updateVars(listVarNames,listVarVals,nbor.name,1)

          #Armamos la parteA
          parteA=updateParteA(nodos,listVarNames,listVarVals,parteA)
          #Agregamos la parteA y B a la matriz
          matrixA.append(parteA)
          matrixB.append(-elementos[nbor.lmntIndex].val)
          #Tras armar la primera ecuación marcamos el supernodo en revisado, para evitar hacer la misma ecuación desde el otro nodo.
          nod.check=True
          getNode(nodos,nod.supernode).check=True
          break

  #Segunda parte de la ecuación de supernodos
  for nod in nodos:
    if nod.supernode!=None:
          #Reiniciamos los valores para la ecuación siguiente:
          for index in range(len(listVarVals)):
            listVarVals[index]=0
          parteA=[]

          #El calculo de la ecuación siguiente se hace en dos fases, suma de las tensiones de nod y nod.super = 0

          #Divisiones de resistencias en supernodo parte izquierda
          lnVcn=nod.name #Nombre variable del nodo actual izquierdo
          lVcn=0 #Variable del nodo actual izquierdo
          #Divisiones de nodos vecinos
          intVal=0
          for lmntIndex in nod.con:
            if elementos[lmntIndex].tpe==2:
              lVcn+=1/elementos[lmntIndex].val

          #Divisiones de resistencias en supernodo parte derecha
          rnod=getNode(nodos,nod.supernode)
          rnVcn=rnod.name #Nombre variable del nodo actual derecho
          rVcn=0 #Variable del nodo actual derecho

          #Divisiones de nodos vecinos
          for lmntIndex in rnod.con:
            if elementos[lmntIndex].tpe==2:
              rVcn+=1/elementos[lmntIndex].val


          for nbor in nod.neighbor:
            Vnn=nod.name
            Vvn=0
            #Para cada nodo
            #Recorremos todos nuestros nodos
            for findN in nodos:
              #Y buscamos el nodo que corresponde a nuestro vecino, no hace falta comprobar que no se opera con el nodo-supernodo porque esta operación omite fuentes
              if nbor.name==findN.name:
                lmntIndex=nbor.lmntIndex
                #Si el nodo vecino no está resuelto es una variable de la ecuación
                if findN.solve!=True:
                  #Si el elemento entre nuestro nodo actual y el nodo vecino es una resistencia
                  if elementos[lmntIndex].tpe==2:
                    Vnn=nbor.name #Nombre de la variable corresponderá al nodo vecino
                    Vvn=-1/elementos[lmntIndex].val #Valor por el que se multiplica nuestro vecino 1/El elemento entre nuestro nodo y el vecino
                    #Actualizamos el valor en la ecuación
                    listVarVals=updateVars(listVarNames,listVarVals,Vnn,Vvn)
                #Si el nodo vecino está resuelto pero el elemento entre estos es una fuente de corriente
                else:
                  #Si el nodo vecino está resuelto, pero el elemento entre este y el nodo actual es una fuente de voltaje
                  intVal+=-findN.volt/elementos[lmntIndex].val
                #Tanto si el nodo vecino está resuelto o no, si tiene un elemento que es fuente de corriente entre estos entonces:
                if elementos[lmntIndex].tpe==4:
                  #Si la corriente indepente apunta anuestro nodo actual la corriente se resta, del caso contrario esta se suma.
                  if elementos[lmntIndex].way==nod.name:
                    intVal+=-elementos[lmntIndex].val
                  else:
                    intVal+=elementos[lmntIndex].val

          #Mismo proceso pero desde el lado derecho, es decir, con los vecinos de nodo con el que se es supernodo o nod.supernode
          for nbor in getNode(nodos,nod.supernode).neighbor:
            Vnn=nod.supernode
            Vvn=0
            #Para cada nodo
            #Recorremos todos nuestros nodos
            for findN in nodos:
              #Y buscamos el nodo que corresponde a nuestro vecino
              if nbor.name==findN.name:
                lmntIndex=nbor.lmntIndex
                #Si el nodo vecino no está resuelto es una variable de la ecuación
                if findN.solve!=True:
                  #Si el elemento entre nuestro nodo actual y el nodo vecino es una resistencia
                  if elementos[lmntIndex].tpe==2:
                    Vnn=nbor.name #Nombre de la variable corresponderá al nodo vecino
                    Vvn=-1/elementos[lmntIndex].val #Valor por el que se multiplica nuestro vecino 1/El elemento entre nuestro nodo y el vecino
                    #Actualizamos el valor en la ecuación
                    listVarVals=updateVars(listVarNames,listVarVals,Vnn,Vvn)
                #Si el nodo vecino está resuelto pero el elemento entre estos es una fuente de corriente
                else:
                  #Si el nodo vecino está resuelto, pero el elemento entre este y el nodo actual es una fuente de voltaje
                  intVal+=-findN.volt/elementos[lmntIndex].val
                #Tanto si el nodo vecino está resuelto o no, si tiene un elemento que es fuente de corriente entre estos entonces:
                if elementos[lmntIndex].tpe==4:
                  #Si la corriente indepente apunta anuestro nodo actual la corriente se resta, del caso contrario esta se suma.
                  if elementos[lmntIndex].way==nod.name:
                    intVal+=-elementos[lmntIndex].val
                  else:
                    intVal+=elementos[lmntIndex].val

          #Actualizamos el valor de la ecuación por el lado derecho y lado izquierdo
          listVarVals=updateVars(listVarNames,listVarVals,lnVcn,lVcn)
          listVarVals=updateVars(listVarNames,listVarVals,rnVcn,rVcn)

          #La ecuación final corresponde a A*Vcn+ [for n in Vnn.lenght(Vnn[n]*Vvn[n])] = -intVal
          #La forma de armar la matriz será en el mismo orden en que se anotaron los elementos de entrada pero omitiendo las variables resueltas
          #For para respetar el orden de ingreso y tener una matriz final coherente
          for N in nodos:
            indiceVnn=0
            for Vn in listVarNames:
              #Si nuestro nodo es el mismo que el de la variable
              if N.name==Vn:
                parteA.append(listVarVals[indiceVnn])
              indiceVnn+=1
          #Al final del proceso agregamos la ecuación en forma matricial por partes
          matrixA.append(parteA)
          #Agregamos la igualdad
          matrixB.append(-intVal)

          #Tras culminar el proceso cerramos el supernodo para no volverlo a operar
          getNode(nodos,nod.supernode).supernode=None
          nod.supernode=None

  return np.array(matrixA),np.array(matrixB),nodos,listVarNames

A,B,listaNodos,listaVariables=getEcx(elementos,nodos)
print(A,B)
for n in listaNodos:
  print("Voltaje en nodo ",n.name," = ",n.volt)

"""#Eliminación Gaussiana con pivote#

## Implementación del método
"""

"""
Forward Elimination
"""
def forwardElimination(A: np.ndarray, b: np.ndarray) -> np.ndarray:
  #n: variable para manejar el tamaño de la matriz
  n = A.shape[0]
  #Se recorren las matrices de forma conveniente
  for i in range(n):
      for j in range(i+1, n):
        #Se compruba que no se encuentre una división por cero en el proceso en cuyo caso la matriz no es invertible; termina la ejecución y retorna las matrices A y b
        if A[i,i]==0:
          print("matriz no invertible")
          return A,b
        #Se calcula el factor cómún por el cual se haga cero en el valor en la posición j de las ecuaciónes inferiores
        factor = A[j, i] / A[i, i]
        # print(factor)
        #Se multiplica este factor por las ecuaciones inferiores para hacer cero ese valor
        for k in range(i, n):
                A[j][k] -= factor * A[i][k]
        b[j] -= factor * b[i]

        # print(A)
        # print(b)
  return A, b

A,B =forwardElimination(A,B)
print(A,B)

#esta función está pensada para matrices tamaño 3, por lo tanto omitirla para otros casos y utilizar el pivoteo a continuación
def backwardSubstitution(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = A.shape[0] #Tamaño de la matriz A
    x = np.zeros(n) #Se crea una matriz de tamaño N llenada con ceros(0)
    #Se resuelven la ecuaciónes de abajo hacia arriba
    #Primero se calcula el primer valor
    x[-1] = b[-1] / A[-1, -1]
    #Se sustituye el valor en las ecuaciones superiores y se resuelve
    for i in range(n-2, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i+1:], x[i+1:])) / A[i, i]
    return x
x=backwardSubstitution(A,B)
print(x[0],x[1],x[2])

# Función de Pivoteo
def pivote(A_aumentada: np.ndarray) -> np.ndarray:
   # Obtiene el número de filas de la matriz A_aumentada
    n = A_aumentada.shape[0]

    # Itera sobre todas las filas de la matriz
    for i in range(n):
        # Encuentra el índice de la fila con el valor absoluto más grande
        # en la columna actual (desde la fila 'i' en adelante)
        fila_maxima = i + np.argmax(np.abs(A_aumentada[i:, i]))

        # Si la fila con el mayor valor absoluto no es la fila actual,
        # intercambia las dos filas
        if fila_maxima != i:
            A_aumentada[[i, fila_maxima]] = A_aumentada[[fila_maxima, i]]

    # Retorna la matriz después de realizar el pivoteo
    return A_aumentada

# Función de Eliminación de Gauss con Pivoteo
def gaussConPivoteo(A: np.ndarray, b: np.ndarray) -> np.ndarray:

    # Crea la matriz aumentada al unir la matriz A y el vector b
    # (b se convierte a una columna antes de unirlo, esto se realiza con el b.reshap(-1,1))
    A_aumentada = np.hstack([A, b.reshape(-1, 1)])

    # Realiza el pivoteo sobre la matriz aumentada
    A_pivoteado = pivote(A_aumentada)

    # Usa forwardElimination para obtener la matriz triangular superior (U) y
    # el vector b actualizado.
    U, b_nuevo = forwardElimination(A_pivoteado[:, :-1], A_pivoteado[:, -1])

    # Una vez que la matriz está en forma triangular superior, se utiliza
    # `backwardSubstitution` para encontrar el vector solución x.
    x = backwardSubstitution(U, b_nuevo)

    # Retorna el vector solución x
    return x

"""## Resultado

Esta celda será utilizada para mostrar el resultado del analisis del circuito, con ambos códigos anexados (nodos y eliminación)
"""

x = gaussConPivoteo(A, B)
print(listaVariables)
print(x)

#Actualizaremos la lista de nodos:
for n in listaNodos:
  ind=0
  for v in listaVariables:
    if n.name==v:
      if n.solve!=True:
        n.volt=x[ind]
    ind+=1

for n in listaNodos:
  print("Voltaje en nodo ",n.name," = ",n.volt)

