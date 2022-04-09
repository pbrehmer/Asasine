## Asasine
Generates audio-visual stream from the solution of the Kuramoto-Sivashinsky partial differential equation. The visual component is a heatmap plot of the solution `U[t, x]`. This same solution at every timestep `U[t, :]` is mapped onto a sine stack, which provides some creative freedom.

The space indices of `U[t, x]` correspond to the sine frequencies of the stack and are (potentially) mapped by some non-linear funtion. Each sine in the stack has an amplitude determined via the actual values of the solution array, which is again put through some mapping; here polynomials `x -> x^7` are one practical option, because they accentuate the non-zero `streamlines` and attenuate the channels inbetween, leading to a more distinct audio spectrum.

### Todo
- figure out: are audio glitches a performance problem or some misunderstanding of `PortAudio.jl`?
  - write up quick function to pre-render audio buffer for some `U[t,x]`
  - smoother interpolation probably doesn't work anymore; think about new one without `last_phases`
- start audio generation and heatmap animation from `U[0, :]`
  - think about renormalizing audio then
- `Makie.jl`:
  - code up own colormap (2 colors for clear L/R channel visual feedback)
  - add left/right spectrum (amp- and frequency-mapped `U[t,x]` slices) to layout
  - add parameter sliders to change parameters on the fly
- investigate divergence of CNAB2 stepping for certain `kx` orderings

### Ideas
