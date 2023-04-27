# Shaft-design-team-6
All user input is within ShaftDesignOptimization.mlx <br />
The user should first verify inputs of shaft specs, material, external loadings, gear specs. <br />
-For shaft geo, the x position that a diameter starts at is the first index then the diameter relation.<br/>
Note that for the first shaft the torque is applied at zero to the gear position and same for the output shaft.<br/>
-For shaft specs, user must input Yielding and Fatigue safety factor, and fillet design factors.<br />
-For material, user inputs Sy, Sut, and rho, these fields can be changed in Shaft.m <br />
-For external loadings, user inputs load multiplier, required torque and rpm. Note that rpm is not used, infinite life was assumed <br />
-For gear specs, user must input the normal pressure angle, normal diametral pitch, helix angle(if spur, = 0), and number of teeth. There are 4 inputs for teeth count and they are sequential from input to output shaft. <br />
Run the .mlx file<br />
The output for fillet radi, final diameter and final weight can be seen first. <br />
User can now look at plots and diameter outputs for all three shafts. Opening plots using the arrow will allow for greater detail and seeing exact points.<br />
For each graph, the x axis represents inches along the shaft. <br />
First the user will see shear stress graphs, followed by a combined bending moment graph.<br />
There is a stress plot for Sigma a and Sigma m <br />
Then, the user can open the diameter figure and see each geometry change in the shaft.<br />
Finally the there is a plot for yielding and fatigue safety facotrs along the shaft. <br />
Process for viewing outputs is the same for all three shafts and the .mlx file will only have to be ran once. <br />
