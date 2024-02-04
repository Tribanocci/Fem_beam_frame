# FEM Frame Calculator
This MATLAB script enables the calculation of the beams stress and displacements due to tensile/compressive and bending loads. The script is user friendly, as the input of the frame geometry can be done using an input.txt file

- [FEM Frame Calculator](#fem-frame-calculator)
  - [Example of Frame Calculation](#example-of-frame-calculation)
    - [Case details](#case-details)
    - [Mathematical formulas](#mathematical-formulas)
    - [Beam data](#beam-data)
    - [Result calculation](#result-calculation)
  - [Use](#use)
    - [Frame geometry](#frame-geometry)
    - [Boundary condition and concetrated loads](#boundary-condition-and-concetrated-loads)


## Example of Frame Calculation

The following frame is going to be analyzed and for each beam the stresses and displacments will be calculated

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.002.jpeg "Καταγραφή")

### Case details

In the above frame we have 8 beams with different directions in the plane. In addition to the nodal forces, we also have distributed forces, while the load due to the same weight has been added to each beam. In addition, the framework was resolved, having as support, bracing at the lower nodes and joints.

### Mathematical formulas
For each beam we have 3 degrees of freedom at each node, so 6 degrees of freedom in total. Therefore the local stiffness matrix k, has dimensions 6x6 with the following form:

![C:\Users\user\AppData\Local\Microsoft\Windows\INetCache\Content.Word\ματ.jpg](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.003.jpeg)

But due to the random direction of each beam, as in lattices, we transform k, so that the projections of the axial and shear forces on the beam enter the global stiffness matrix. Therefore displacements are calculated in global coordinates.

Distributed and concentrated forces are entered according to global coordinates. Furthermore, the distributed forces are transformed into epicondials through the integral.

![C:\Users\kfyti\Downloads\A3\Pic1.png](doc\Pic1.png)

Although in the elements between nodes (4,5) and (5,6), they receive a distributed force not perpendicular to the neutral axis of the beam, no further transformation is needed. The reason this happens is because the distributed forces are transformed into epinodal forces in the global system and so we can proceed with the solution for the displacements of the frame, which are done in the global coordinate system.

The loading due to the same weight is included with the corresponding epicondicular forces and moments which produces a constant distributed force along the length of the beam. But if we have a vertical beam, the loading of the weight consists of two epinodal forces in the y direction, divided between the two nodes.

*The above calculations are done with the following piece of code:*

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.004.jpeg)

### Beam data

The beam used for the frame is an i-beam, corresponding to the identification IPE-200. The cross-section is shown in the figure below and has the characteristics shown in the next table.



![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.005.jpeg "ibeam")

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.006.jpeg "ipe200")

From these data we keep the area, the moment of inertia about the strong axis (for bending) and the density.

Area=28.5 cm2   I=1943 cm4   ρ=22,4 kg/m

While the modulus of elasticity was considered equal to Ε=200 GPa

From the solution of the framework, the displacements of the framework are obtained, where if K is the global stiffness matrix, and F the epinodal forces in the global system, the displacements in global coordinates are obtained as Dglobal=K-1F

ΓTo find the forces at each node in the local system, the linear transformation is done 

f=k\*T\*Dglobal

Also, if in a beam we have the presence of a distributed force, then the equivalent epinodal forces must be removed from the above calculation (which are also found in the global system). This is done by the following act:



f=k\*T\*Dglobal-Τ\*Fdist

*The above calculations are done with the following piece of code:*

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.007.jpeg "Φ")

### Result calculation

The shear force and bending moment diagrams were obtained from finite element theory for the beam, where from the cubic shape function N and from the relations

![C:\Users\kfyti\Downloads\A3\Pic1.png](doc\pic11.png)

The axial force diagrams were obtained from the local forces of each element in the x direction. In the beams where we have no distributed force in the x-axis, the construction of the diagram is simple, since we have constant axial stress within the beam. But for the two beams at the top of the frame, the distributed force has an x-axis component as well. Therefore, the diagram of the axial forces has a linear form and the values ​​of the local forces in the x direction, give the two extreme values.



*The above calculations are done with the following piece of code:*

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.008.jpeg "πλ")

Summarizing, the code [*beam\_frame.m*] having as input the coordinates of the nodes, the epinodal and distributed forces on each edge (beam) and the characteristics of the cross section, gives the displacements of the nodes, the local forces on each node and the NQM diagrams for each beam.









1<sup>η</sup> Δοκός

![C:\Users\user\AppData\Local\Microsoft\Windows\INetCache\Content.Word\1a.jpg](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.009.jpeg)






2<sup>η</sup> Δοκός

![C:\Users\user\Desktop\2a.JPG](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.010.jpeg)








3<sup>η</sup> Δοκός



![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.011.jpeg "3a")


4<sup>η</sup> Δοκός



![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.012.jpeg "4a")

5<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.013.jpeg "5a")



6<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.014.jpeg "6a")





7<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.015.jpeg "7a")




8<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.016.jpeg "8a")



Τα ακόλουθα διαγράμματα για κάθε δοκό είναι για πάκτωση 

1<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.017.jpeg "1b")



2η Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.018.jpeg "2b")



3<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.019.jpeg "3b")



4<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.020.jpeg "4b")




5<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.021.jpeg "5b")




6<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.022.jpeg "6b")




7<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.023.jpeg "7b")



8<sup>η</sup> Δοκός

![](doc\Aspose.Words.53f010bd-91cb-4803-ad61-28c168fb6218.024.jpeg "8b")




From the above analysis we conclude that when we have joint support the vertical beams at the base have a greater stress due to bending moments.

The maximum bending moment for joint support appears in the beam of the right base at node 8 with a value of M_8=101.16 kNm, while the corresponding bending moment for compaction is M_8=75.75 kNm.

Also, the same picture applies to the roof beams, with those corresponding to the frame with the joint showing higher bending moments than those of the compaction. However, the maximum bending flow occurs in both cases at the base of the frame.

## Use

In order to input the data:

### Frame geometry

In the *frame.txt*, in each line pass 4 numbers, with the following order seperated by space:
node_start_index node_final_index Length_in_X_direction Length_in_Y_direction

### Boundary condition and concetrated loads

In the *frame_dist_force.txt*, in each line pass 4 numbers, with the following order seperated by space:
node_index, Force_in_X_diection, Force_in_Y_diection, Moment_in_node, Fixed_in_X_direction, Fixed_in_Y_direction, Fixed_in_moment


node_start_index node_final_index Force_start_value Force_final_value