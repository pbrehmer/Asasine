# Asasine

Generates audio-visual stream from the solution of the Kuramoto-Sivashinsky partial differential equation. The visual component is a heatmap plot of the solution `U[t, x]`. This same solution at every timestep `U[t, :]` is mapped onto a sine stack, which provides some creative freedom.

The space indices of `U[t, x]` correspond to the sine frequencies of the stack and are (potentially) mapped by some non-linear funtion. Each sine in the stack has an amplitude determined via the actual values of the solution array, which is again put through some mapping; here polynomials `x -> x^7` are one practical option, because they accentuate the non-zero `streamlines` and attenuate the channels inbetween, leading to a more distinct audio spectrum.

## Todo

- add interactive elements
  - fix button: fixes `U[x,t]` slice and halts evolution
  - restart button
- generate audio buffer via IFFT for improved performance
- add export scripts: export audio/video/image separately
- new colormap? (2 colors for clear L/R channel visual feedback)
- stop divergence of CNAB2 stepping for certain `kx` orderings
