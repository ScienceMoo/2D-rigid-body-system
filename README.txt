# Rigid Body Collision and Contact

Steps and Objectives

1. Rotational Rigid Body Dynamics
- This was fairly easy to do and seems to work if LCP solve is off.
- When using LCP the torque from object collision is calculated differently and may be wrong, I am not sure

2. Hierarchical Collision Detection
- I got this to work successfully by simply checking both child nodes one after another. It is definitely much faster.

3. LCP solve
- I am still not sure if this is working correctly, even after implementing all the optimizations in the next steps
- It works with no bugs only for the simple test examples

4. Optimization
- Lambda seemed to converge quickly enough for small systems, but I was not able to test my code for larger systems
- Randomizing the order of the inner loop improved convergence of lambda
- I implemented randomization using a shuffled list of indeces

5. Warm Starts
- I implemented this using hashmaps

6. Demonstration movies, Novel Scene with Interesting Movie
- Unfortunately I could not show certain scenes due to bugs that I couldn't figure out, but I did make a novel demonstration movie and showed a few scenes with and without LCP solve. 
- Furthermore, I was not able to record and simulate at the same time while doing LCP due to lack of memory space, I must have done a bad job allocating memory... so I recorded the scenes using the iOS screen recorder which is still pretty high quality. I edited them with iMovie.
- I wish I could have spent more time on this assignment to figure everything out.