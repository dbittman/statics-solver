# statics-solver
Solver for statics problems, given boundary conditions, using the Poisson equation.

This projects the work for my physics senior thesis at UC Santa Cruz. I plan to implement varying
ways to actually solve the resulting linear algebra problem after discretizing a region for the Poisson equation.
I've planned to do Jacobi iterations, successive-over-relaxation, and possibly a Fourier transform method using the FFT.
Then I'll compare performance. I'm also going to make it multithreaded and use SIMD.

Fun stuff! It actually kinda works right now, though it's slow and the region conditions are hard-coded. But you can see it work:
``make && ./main 600 && python grad.py``.
