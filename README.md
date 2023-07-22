# Asasine

Generates audio and visual data from the solution of the Kuramoto-Sivashinsky (KS) partial differential equation

$$\frac{\partial u}{\partial t} + \frac{\partial^2 u}{\partial x^2} + \frac{\partial^4 u}{\partial x^4} + \frac{1}{2} \Big(\frac{\partial u}{\partial x}\Big)^2 = 0 .$$

The visual component is a heatmap plot of the solution function $u(t, x)$, while the audio component is the solution mapped onto a frequency spectrum that is written to an audio buffer at each time $t$. This audio-visual feed can be streamed in real time or rendered to a file.

The sonification process is implemented as follows: At each time step, from the KS solution vector `u` a collection of values is picked via an index vector `freqidx` with corresponding frequencies `freqs`. The solution values determine the sine wave amplitudes at the respective frequencies. Here, much of the sound characteristics and creative freedom resides in the choice of `freqs`. Additionally, we apply a mapping `ampmap` that further modifies the amplitudes of the produced sine waves; we often times use monomials of odd order, e.g., `ampmap = a -> a^5` to sharpen the features of the produced frequency spectrum and still retain negative values ($u(t, x)$ takes on both negative and positive real values). In simplified (code) terms, we generate the audio buffer at each time step via:

```julia
buf = zeros(buflen)
for (i, ν) in zip(freqidx, freqs)
    buf .+= ampmap(u[freqidx[i]]) .* sin.(2pi * buflen * ν)
end
```

In order to create an interesting stereo image, the KS solution is translated to a stereo buffer by mapping all negative and positive values of `u` to the left and right channel, respectively.

For more details, check the notes in the [introductory notebook](/examples/introduction.ipynb) as well as the provided [examples](/examples/).

## Sample outputs

![L128_1](https://github.com/pbrehmer/Asasine/assets/62562093/09599fc6-6796-4c74-b1d0-3811c857c73c)

https://github.com/pbrehmer/Asasine/assets/62562093/12be5335-50ec-4fb0-9d9c-f89316e291ab

![L256_1](https://github.com/pbrehmer/Asasine/assets/62562093/21b960f9-f210-4dc1-8229-11d72f8af245)

https://github.com/pbrehmer/Asasine/assets/62562093/708c931b-094a-43c5-b908-2e7d0a1b55d6

![L128_2](https://github.com/pbrehmer/Asasine/assets/62562093/57fc027b-a93d-4347-85b9-ce4f1a47b704)

https://github.com/pbrehmer/Asasine/assets/62562093/fccf0a0d-4552-45a1-af34-3017ef682f7c

https://github.com/pbrehmer/Asasine/assets/62562093/8ab0b70b-41a6-4a5b-971f-e21d810e4236

![L32_1](https://github.com/pbrehmer/Asasine/assets/62562093/84388bd7-ad26-4f17-96c9-7602d2187747)

https://github.com/pbrehmer/Asasine/assets/62562093/fbcf4d37-eda0-443f-93db-93f712121480

## References

The idea for this project was sparked by a blog post series by John Carlos Baez [1] on the Kuramoto-Sivashinsky equation. The CNAB2 Julia implementation was largely inspired by Mathab Lak's code [2], which is also used in Ref. [1].

[1] John Carlos Baez, [The Kuramoto–Sivashinsky Equation](https://johncarlosbaez.wordpress.com/2021/10/17/conjectures-on-the-kuramoto-sevashinsky-equation/)\
[2] Mathab Lak, [Test case for PDEs: Kuramoto-Sivashinsky](https://online.kitp.ucsb.edu/online/transturb17/gibson/html/5-kuramoto-sivashinksy.html)

## Todo

- add interactive elements
  - fix button: fixes `U[x,t]` slice and halts evolution
  - restart button
- generate audio buffer via IFFT for improved performance
- investigate divergence of CNAB2 stepping for certain wave number orderings
