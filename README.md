# Asasine
Generates audio-visual stream from the solution of the Kuramoto-Sivashinsky partial differential equation. The visual component is a heatmap plot of the solution `U[t, x]`. This same solution at every timestep `U[t, :]` is mapped onto a sine stack, which provides some creative freedom.

The space indices of `U[t, x]` correspond to the sine frequencies of the stack and are (potentially) mapped by some non-linear funtion. Each sine in the stack has an amplitude determined via the actual values of the solution array, which is again put through some mapping; here polynomials `x -> x^7` are one practical option, because they accentuate the non-zero `streamlines` and attenuate the channels inbetween, leading to a more distinct audio spectrum.

# Todo
- circumvent `on(...)` multiple execution in while-loop
- add interactive elements
  - start/step/stop slider for frequency indices
  - restart button
- maybe generate audio buffer via IFFT for improved performance
- add record function
  - export audio/video separately or fuse into one `.mp4` file?
- start audio generation and heatmap animation from `U[0, :]`
  - think about renormalizing audio then
- code up own colormap (2 colors for clear L/R channel visual feedback)
- investigate divergence of CNAB2 stepping for certain `kx` orderings
