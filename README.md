# Asasine

Generates audio and visual data from the solution of the Kuramoto-Sivashinsky (KS) partial differential equation

$$u_t + u_{xx} + u_{xxxx} + \frac{1}{2} u_x^2 = 0 .$$

The visual component is a heatmap plot of the solution function $U(t, x)$, while the audio component is the solution mapped onto a frequency spectrum that is written to an audio buffer at each time $t$. This audio-visual feed can be streamed in real time or rendered to a file.

The sonification process is implemented as follows: At each time step, from the KS solution vector `u` a collection of values is picked via a vector frequency indices `freqidx` that corresponds to the frequencies `freqs`. Additionally, we apply a mapping `ampmap` modifies the amplitudes of the produced sine waves; we often times use monomials of odd order, e.g., `ampmap = a -> a^5` to sharpen the features of the produced frequency spectrum and still retain negative values ($U(t, x)$ takes on both negative and positive real values). Altogether, we in principle generate an audio buffer `buf` of length `buflen` at each time step via:

```julia
for (i, ν) in zip(freqidx, freqs)
    buf .+= ampmap(u[freqidx[i]]) .* sin.(2pi * buflen * ν)
end
```

In order to create an interesting stereo image, the KS solution is translated to a stereo buffer by mapping all negative and positive values of `u` to the left and right channel, respectively.

For more details, check the notes in the [introductory notebook](/examples/introduction.ipynb) as well as the provided [examples](/examples/).

## Todo

- add interactive elements
  - fix button: fixes `U[x,t]` slice and halts evolution
  - restart button
- generate audio buffer via IFFT for improved performance
- add export scripts: export audio/video/image separately
- investigate divergence of CNAB2 stepping for certain `kx` orderings
