Read me
=======

  output         	 = output file
  L              	 = system size
  s,a,epsilon     	 = reaction terms coefficients
  diffprex,diffpak       = diffusion coefficients
  dt             	 = time step
  totaltime      	 = duration of the run (after equilibration)
  initial_time   	 = set t=0
  interval_record 	 = time interval between two successive recordings
  dx              	 = space discretization

      L       s          a     epsilon    diffprex    diffpak       dt       totaltime   initial_time   interval_record    seed        dx

RUN  1e0      0.03       1.5     0.07       0.001      0.5       1e-4        1e3          0             6e0               $RANDOM      2e-2



 
To compile

gcc -O3 Finite_difference_simulation.c -o RUN -lm


