# Gaussian Elimination to solve a circuit by nodal method
# Circuit modelation and Gauss Elimination for solve circuits by node method. Dec 2023.

"""
I expanded the scope of a work designed for numerical methods. The code didn't was optimize but it could be use as a reference, 
the reach of this code was proved with examples circuitos, could work in cases of resisntors, voltage sources and current sources 
(as long as they are independent) and supernode cases. Maybe it could not work properly in cases with supernodes with 3 or more nodes,
and could get difference with real results. Anyway it could be use as reference. There is the file .ipynb ideal for work in collab, also 
there is the file .py with the code and comments, both with spanish descriptions by now.

By
Leandro Maureira.


# Analysis and get nodes
"""

import numpy as np
import json

#Class Definition
class interlmnt:
    def __init__(self,neighbor,lmntIndex):
      self.name=neighbor #Neighbor node name
      self.lmntIndex=lmntIndex #index of an element in between of two nodes

class node:
    def __init__(self,name,con,ref):
      self.name=name
      self.con=con #Array with all elements conected to a node
      self.ref=ref #Boolean, true if it's and reference node
      self.way=[] #Path of current. #+1 get in #-1 get out #[Wasn't used]
      self.volt=None #Voltage of the node, None if isn't defined/solved
      self.solve=False #True:node solved False:node unsolved
      self.neighbor=[] #Neighbor list in tuples of interlmnt class
      self.supernode=None #If this node it's a supernode it says wich with node it is
      self.check=False #False by default #In case of supernode it say if this node was checked

class lmnt:
    def __init__(self,ind,tpe,val,way):
      self.ind=ind #Indice o numero del elemento #Index/Name of an element

      #Values for tpe(Type of element)
      #0 Dependent voltage source #Not implemented
      #1 Independent voltage source #Working
      #2 Resistors #Working
      #3 Dependent current source #Not implemented
      #4 Independent current source #Working
      self.tpe=tpe

      self.val=val #Value

      self.way=way #To wich node goes the current/To wich node aim the negative signal

#JSON example to model a circuit
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



data=json.loads(exampleD)

#Read JSON and save data en the lists
#Return two lists, first one with elements(sources, resistors,etc) and second one the nodes
def modCirc(jsonfile):
  elementList=[]
  nodeList=[]
  for i in data["lmnts"]:
    element = lmnt(i[0],i[1],i[2],i[3])
    elementList.append(element)
    #print(elemento.val)
  for j in data["nodes"]:
    nod = node(j["name"],j["con"],j["ref"])
    nodeList.append(nod)

  return elementList,nodeList

#Simulation of circuit
elements,nodes=modCirc(data)

#Helping functions
#Update the value of the multiplication that goes with the variable
def updateVars(listVarName,listVarVal,VarName,VarVal):
  index=0
  for name in listVarName:
    if name==VarName:
      listVarVal[index]+=VarVal
    index+=1
  return listVarVal

#Check if a value it's already in listVarNames
def varDefPre(listVarNames,name):
  for n in listVarNames:
    if n==name:
      return True
      break
  return False

#Check if an specific node it's a reference node
def isRefNode(name,nodeList):
  for n in nodeList:
    if n.name==name and n.ref==True:
      return True
  return False


#It form the partA of the matrix in 'First In' order. Recieves node list, variable names list, values list and current part A
def updatePartA(nodes,listVarNames,listVarVals,partA):
  for N in nodes:
        indiceVnn=0
        for Vn in listVarNames:
          #If our node it's the same that the variable
          if N.name==Vn:
            partA.append(listVarVals[indiceVnn])
          indiceVnn+=1
  return partA

#Look the 'looking' node and return that node
def getNode(nodes,looking):
  for n in nodes:
    if looking==n.name:
      return n
  return None

#Get ecuations of voltage from the circuit and return in two matrixs A,B
#Return Matrix A,Matrix B,node list updated with values of solved nodes(Reference node and sources to earth), and return variable list in order
#Working with voltage sources, current sources, resistors, doesn't proved with supernodes with more than two nodes
def getEcx(elements,nodes):
  #Nodal method, going throw each node

  #Look for neighbor nodes
  for CurrentN in nodes:
    #We solved reference node
    if CurrentN.ref==True:
      CurrentN.volt=0
      CurrentN.solve=True
    for N in nodes:
      #Check if isn't comparing a node with itself
      if CurrentN.name!=N.name:
        #Going throw conections of each node
        for currentCon in CurrentN.con:
          for Con in N.con:
            #If it has a common conection it's a neighbor node
            if currentCon==Con:
              #We save the neighbor node with the element in between
              CurrentN.neighbor.append(interlmnt(N.name,currentCon))

  listVarNames=[] #List with names of variable to calculate
  listVarVals=[] #List with matricial values of each variable

  for nod in nodes:
    #For not reference nodes
    if nod.ref!=True:
      #Check if there's voltage sources in this node
      #Going throw each neighbor
      for nb in nod.neighbor:
        #If there is voltage sources this node value it's equal to the value of the voltage just in case it's connected to earth
        if elements[nb.lmntIndex].tpe==1 and isRefNode(nb.name,nodes)==True:
          #Consider where goes the negative of the source
          if elements[nb.lmntIndex].way==nb.name:
            nod.volt=elements[nb.lmntIndex].val
          else:
            nod.volt=-elements[nb.lmntIndex].val
          #Tick as solved
          nod.solve=True
          #If we put this variable in the list we should take off in this case
          listVarNames.remove(nod.name)
         #Now it is solved we break the loop cuase we have no more to do in this node.
          break
        #If the source of voltage isn't connected to earth then we have a supernode with the node that we have this source in common
        #Check if this node isn't solved
        elif elements[nb.lmntIndex].tpe==1 and isRefNode(nb.name,nodes)==False and nod.solve!=True:
          nod.supernode=nb.name
        #If the variable have no direct solution so there is an incongnit and now it's a variable from the equation
        #Check that the variable wasn't add
        if varDefPre(listVarNames,nod.name)==False and nod.solve!=True:
          listVarNames.append(nod.name)
          listVarVals.append(0)
  print(listVarNames)

  #Restart the loop with solved variables and make the matrixs
  matrixA=[]
  matrixB=[]
  #Analysis, node by node without solved nodes and supernodes
  for nod in nodes:
    #Re-initialize values of unknowed variables for the current equation in this node
    for index in range(len(listVarVals)):
      listVarVals[index]=0

    if nod.solve!=True and nod.supernode==None:
      partA=[]
      #In nodal method if there's no voltage sources the outputs currents are equal to 0
      #We make the equation with neighbor nodes
      #Resistors divisions
      nVcn=nod.name #Name of the variable in the current node
      Vcn=0 #Value of the variable of the current node
      #Divisions in neighbor nodes
      intVal=0
      for lmntIndex in nod.con:
        if elements[lmntIndex].tpe==2:
          Vcn+=1/elements[lmntIndex].val

      #Look for the neighbor node except supernodes
      if nod.supernode==None:
        for nbor in nod.neighbor:
          Vnn=nod.name
          Vvn=0
          #We going throw each node by each node
          for findN in nodes:
            #For each node we define a list of used elements (?) [Deleted line]
            #We look for the node wich it's our neighbor
            if nbor.name==findN.name:
              lmntIndex=nbor.lmntIndex
              #If this node wasn't solve it's an incognit of the equation
              if findN.solve!=True:
                #If the element in between our current node and the neighbor node it's a resistor
                if elements[lmntIndex].tpe==2:
                  Vnn=nbor.name #The variable incognit it's the neighbor node
                  Vvn=-1/elements[lmntIndex].val #The value of the multiplication is 1/value of lmnt in between
                  #Update the value in the equation
                  listVarVals=updateVars(listVarNames,listVarVals,Vnn,Vvn)
              #If the neighbor node it's unsolved but the lmnt in between it's a current source
              else:
                #If the neighbor node it's solved but the lmnt in between it's a voltage source
                intVal+=-findN.volt/elements[lmntIndex].val
              #Tanto si el nodo vecino está resuelto o no, si tiene un elemento que es fuente de corriente entre estos entonces:
              #Wheter the node it's solved or not, if there's a current source then:
              if elements[lmntIndex].tpe==4:
                #Si la corriente indepente apunta anuestro nodo actual la corriente se resta, del caso contrario esta se suma.
                #If the independent current source aim to our current node the current it's substracted , opposite case it's added.
                if elements[lmntIndex].way==nod.name:
                  intVal+=-elements[lmntIndex].val
                else:
                  intVal+=elements[lmntIndex].val

      #Update the value in equation
      listVarVals=updateVars(listVarNames,listVarVals,nVcn,Vcn)

      #The final equation it's A*Vcn+ [for n in Vnn.lenght(Vnn[n]*Vvn[n])] = -intVal
      #We will make the matrix in the same order as was read but without the solved variables
      #That's why we use 'for'.
      for N in nodes:
        indiceVnn=0
        for Vn in listVarNames:
          #If neighbor node is the same than the variable
          if N.name==Vn:
            partA.append(listVarVals[indiceVnn])
          indiceVnn+=1
      #Finally we put the equation in matrix form by parts
      matrixA.append(partA)
      #Put the equality
      matrixB.append(-intVal)

  #Supernodes exceptions; a source of voltage wich it's not connected to earth
  #First part
  #In this case there will be two equations; the substract a node from the other it's equal to the voltage in between
  for nod in nodes:
    #If the node is part of a supernode and wasn't check before
    if nod.supernode!=None and nod.check==False:
      #Restart our values in the incognit variables for the equation in the current node
      for index in range(len(listVarVals)):
        listVarVals[index]=0
      partA=[]
      #Look for the neighbor node wich is part of this supernode
      for nbor in nod.neighbor:
        #Check if it's the correct element in between(voltage source)
        if nod.supernode==nbor.name and elements[nbor.lmntIndex].tpe==1:
          #Then the supernode it's nod + nod.supernode

          #Define the first equation, nod-nbor=voltage of the middle element, this means (1 -1) (Vx), negative for the node who isn't aimed by negative in the source.
          #Opposite case (-1 1) (Vx)
          if elements[nbor.lmntIndex].way==nod.name:
            listVarVals=updateVars(listVarNames,listVarVals,nod.name,1)
            listVarVals=updateVars(listVarNames,listVarVals,nbor.name,-1)
          #Opposite case
          else:
            listVarVals=updateVars(listVarNames,listVarVals,nod.name,-1)
            listVarVals=updateVars(listVarNames,listVarVals,nbor.name,1)

          #Make partA matrix
          partA=updatePartA(nodes,listVarNames,listVarVals,partA)
          #Add partA and B
          matrixA.append(partA)
          matrixB.append(-elements[nbor.lmntIndex].val)
          #Tick reviwed (To avoid check from the other side node )
          nod.check=True
          getNode(nodes,nod.supernode).check=True
          break

  #Second part in the supernode equations
  for nod in nodes:
    if nod.supernode!=None:
          #Restart values for the next equation
          for index in range(len(listVarVals)):
            listVarVals[index]=0
          partA=[]

          #The equation will be make in two steps, voltage sums and nod.super = 0

          #Divisions of resistors on the left side of supernode
          lnVcn=nod.name #Variable name of current left node
          lVcn=0 #Variable name of current left node

          #Divisions of neighbor nodes
          intVal=0
          for lmntIndex in nod.con:
            if elements[lmntIndex].tpe==2:
              lVcn+=1/elements[lmntIndex].val

          #Divisions of resistors on the right side of supernode
          rnod=getNode(nodes,nod.supernode)
          rnVcn=rnod.name #Variable name of current right node
          rVcn=0 #Variable name of current right node

          #Divisions of neighbor nodes
          for lmntIndex in rnod.con:
            if elements[lmntIndex].tpe==2:
              rVcn+=1/elements[lmntIndex].val


          for nbor in nod.neighbor:
            Vnn=nod.name
            Vvn=0
            #For eachnode we go throw each node
            for findN in nodes:
              #Look for the node wich it's the neighbor, this part avoid source then we will not check that
              if nbor.name==findN.name:
                lmntIndex=nbor.lmntIndex
                #If the neighbor isn't solved is a variable in the equation
                if findN.solve!=True:
                  #If the element in the middle is a resistor
                  if elements[lmntIndex].tpe==2:
                    Vnn=nbor.name #Name of varaible in neighbor node
                    Vvn=-1/elements[lmntIndex].val #Value in the multiplication 1/element in between
                    #Update the value in the equation
                    listVarVals=updateVars(listVarNames,listVarVals,Vnn,Vvn)
                #If the neighbor node is solved but the element in between is a voltage source
                else:
                  #If the node is solved but the element in the middle is a voltage source
                  intVal+=-findN.volt/elements[lmntIndex].val
                #As if the neighbor node is solved as is not, if there is a current source in between then:
                if elements[lmntIndex].tpe==4:
                  #If the independent current source aim to the current node this is substracted, opposite case there is added
                  if elements[lmntIndex].way==nod.name:
                    intVal+=-elements[lmntIndex].val
                  else:
                    intVal+=elements[lmntIndex].val

          #Same process but in the right side; nod.supernode
          for nbor in getNode(nodes,nod.supernode).neighbor:
            Vnn=nod.supernode
            Vvn=0
            #For eachnode we go throw each node
            for findN in nodes:
              #Look for the node wich it's the neighbor
              if nbor.name==findN.name:
                lmntIndex=nbor.lmntIndex
                #If the neighbor isn't solved is a variable in the equation
                if findN.solve!=True:
                  #If the element in the middle is a resistor
                  if elements[lmntIndex].tpe==2:
                    Vnn=nbor.name #Name of varaible is neighbor node
                    Vvn=-1/elements[lmntIndex].val #Value in the multiplication 1/element in between
                    #Update the value in the equation
                    listVarVals=updateVars(listVarNames,listVarVals,Vnn,Vvn)
                #If the neighbor node is solved but the element in between is a voltage source
                else:
                  #If the node is solved but the element in the middle is a voltage source
                  intVal+=-findN.volt/elements[lmntIndex].val
                #As if the neighbor node is solved as is not, if there is a current source in between then:
                if elements[lmntIndex].tpe==4:
                  #If the independent current source aim to the current node this is substracted, opposite case there is added
                  if elements[lmntIndex].way==nod.name:
                    intVal+=-elements[lmntIndex].val
                  else:
                    intVal+=elements[lmntIndex].val

          #Updated the value in the equation, right side and left side
          listVarVals=updateVars(listVarNames,listVarVals,lnVcn,lVcn)
          listVarVals=updateVars(listVarNames,listVarVals,rnVcn,rVcn)

          #Final equation will be A*Vcn+ [for n in Vnn.lenght(Vnn[n]*Vvn[n])] = -intVal
          #We will make the matrix in the same order as was read but without the solved variables
          #That's why we use 'for'.
          for N in nodes:
            indiceVnn=0
            for Vn in listVarNames:
              #If neighbor node is the same than the variable
              if N.name==Vn:
                partA.append(listVarVals[indiceVnn])
              indiceVnn+=1
          #Finally we put the equation in matrix form by parts
          matrixA.append(partA)
          #Put the equality
          matrixB.append(-intVal)

          #Then close the node for no operate again
          getNode(nodes,nod.supernode).supernode=None
          nod.supernode=None

  return np.array(matrixA),np.array(matrixB),nodes,listVarNames

A,B,nodeList,variableList=getEcx(elements,nodes)
print(A,B)
for n in nodeList:
  print("Voltage in node ",n.name," = ",n.volt)

"""#Gaussian Elimination

## Method
"""

"""
Forward Elimination - traditional implementation
"""
def forwardElimination(A: np.ndarray, b: np.ndarray) -> np.ndarray:
  #n: #Size of matrix
  n = A.shape[0]
  #going conveniently in the matrix
  for i in range(n):
      for j in range(i+1, n):
        #Check isn't a division by cero in the process, then isn't invertible
        if A[i,i]==0:
          print("matrix isn't invertible")
          return A,b
        #Calculate the common factor to make cero the value in position j of below equations
        factor = A[j, i] / A[i, i]
        # print(factor)
        #Multiplicate this factor for the belows equations to make 0
        for k in range(i, n):
                A[j][k] -= factor * A[i][k]
        b[j] -= factor * b[i]

        # print(A)
        # print(b)
  return A, b

A,B =forwardElimination(A,B)
print(A,B)

#Just for size 3 matrix, avoid in other cases and use pivot method
def backwardSubstitution(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = A.shape[0] #Size of matrix A
    x = np.zeros(n) #fill with ceros the matrix x size n
    #Solve from bottom to top
    #First value
    x[-1] = b[-1] / A[-1, -1]
    #Substitution on upper equations
    for i in range(n-2, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i+1:], x[i+1:])) / A[i, i]
    return x
x=backwardSubstitution(A,B)
print(x[0],x[1],x[2])

# Pivot method
def pivot(A_increased: np.ndarray) -> np.ndarray:
   # Get rows count in A_increased
    n = A_increased.shape[0]

    # Loop in all rows
    for i in range(n):
        # Find index of row with bigger absolut value
        # In current column (from i)
        max_row = i + np.argmax(np.abs(A_increased[i:, i]))

        # If the row with bigger absolut value isn't the current
        # Swap the two rows
        if max_row != i:
            A_increased[[i, max_row]] = A_increased[[max_row, i]]

    # Return matrix after do the pivote
    return A_increased

# Gauss with pivote
def pivotedGauss(A: np.ndarray, b: np.ndarray) -> np.ndarray:

    # Define the increased matrix A+b
    # b will be a collumn before add, so use b.reshap
    A_increased = np.hstack([A, b.reshape(-1, 1)])

    # Pivote on increased matrix
    A_pivoted = pivot(A_increased)

    # Forward Elimination to get upper triangular (U) and updated b
    U, b_nuevo = forwardElimination(A_pivoted[:, :-1], A_pivoted[:, -1])

    # Then use backward substitution to find the solution vector x
    x = backwardSubstitution(U, b_nuevo)

    # Return x solution
    return x

"""## Results
"""

x = pivotedGauss(A, B)
print(variableList)
print(x)

#Update node list
for n in nodeList:
  ind=0
  for v in variableList:
    if n.name==v:
      if n.solve!=True:
        n.volt=x[ind]
    ind+=1

for n in nodeList:
  print("Voltage in nodo ",n.name," = ",n.volt)
