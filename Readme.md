redist
====
Overview
---
redist is a small utility program that modifies the velocity output generated by [NPSAT](https://gwt.ucdavis.edu/research-tools-and-applications/npsat-engine) so that it can be used by [Ichnos](https://github.com/giorgk/ichnos).

When the NPSAT is executed in a multi-core mode each core prints the locally owned velocity field in its own file. NPSAT is based on [deal.II](https://www.dealii.org/) which in turn uses [p4est](https://www.p4est.org/) to deal with the distribution of degrees of freedom (dofs) to different processors. This results in a distribution of dofs which look like the following:

![Typical distribution of degrees of freedom from p4est](Example_data/dofs_example1_v2.png)

The coloured dots correspond to nodes where the velocity is known from the finite element solution. The dots of the same color belong to the same processor. We can see that p4est distribution creates subdomains that are very irregular. It is also very likely that one domain may contaisn dofs that ther are totally disconnected in the physical space. While this is acceptable when solving linear systems, it causes quite a few problems on particle tracking codes.

Therefore the purpose of this utility program is to redistrube the dofs into a regular pattern that can be used by [Ichnos](https://github.com/giorgk/ichnos). In the above example the pattern is to split the nodes into 4 ractangular domains which are shown here with blue lines. These are called __Actual__ domains in the redist terminology. In particle tracking, when a particle exits the subdomain it is very usefull to know in which subdomain should go. We have found also that it is very usefull each processor to know the velocities for a layer of points around the subdomain. Therefore redist expects a number of __Expanded__ polygons which contain the velocity points to be printed into file.       

How to use
---
First prepare an input file with the following format
```
Nproc prefix leadingZeros suffix useGraph iter
Ndata NprintData Col1 Col2 ...
Ndigit1 Ndigit2 Ndigit3 ...
ActualDomain ExpandedDomain
NewPrefix 
```
where 

`Nproc` is the number of domains in the FEM simulation </br>
`prefix` is the prefix of the velocity file names </br>
`leadingZeros` is the number of zero padding in the velocity file name </br>
`suffix` is the file extension </br>
`useGraph` if this is 1 will read the graph files</br>
`iter` THis is the number of iterations for the buffer zone. A value of 2 is a good choice. It is used only if `useGraph = 1`</br> 
`Ndata` is the number of data in the velocity files without the `x,y,z` coordinates </br>
`NprintData` is the number of data to print. This number equal or less than ` Ndata` </br>
`Col1 Col2 ...` are the column indices of the Ndata in the order that should be printed </br>
`Ndigit1 Ndigit2 Ndigit3 ...` is the number of decimals in the new file</br>
`ActualDomain` The filename that containts the polygons of the actual domains</br>
`ExpandedDomain` The filename that containts the polygons of the expanded domains</br>
`NewPrefix` The new prefix. The filename will have the same leading zeros and extension as the initial files

For example</br>
```
4 example1_ 4 .vel
6 5 5 6 1 2 3 
2 2 6 6 6
example1_actual_dom.dat example1_expanded_dom.dat
examp1e1_new_
```
or
```
6 output2/TuleRiver2_ 4 .vel 1 2
6 5 5 6 1 2 3 
2 2 6 6 6
Tule_actual_dom.dat Tule_extended_dom.dat
output2/TuleRiver2new_
```

For the first input file the FEM simulation was carried out on 4 processors.</br>
The second line indicates that the initial files contain 6 columns of data but only 5 of them will be printed. The columns 5 and 6 will be printed with 2 decimal digits and the columns 1 2 3 with 6 decimal digits. 


In the figure above the left lower corner has coordinates (-5000, -2500) and the upper right corner (5000, 2500). The actual and expanded domain files that correspond to the above figure with blue and red rectangles respectively can be found under this repository

Then run the code using as many cores as the number of polygons in the actual and expanded domain files. For this example this is 4 i.e.
```
cd Example_data
C:\"Program Files\Microsoft MPI"\Bin\mpiexec.exe -n 4 ..\Debug\redist.exe example1_input.txt
```