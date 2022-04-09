using PortAudio
using AuditoryStimuli, Unitful, Plots, Pipe, DSP

# start Julia REPL
test = 1
# set up sink
# devices = PortAudio.devices()
# println(devices)
sample_rate = Float64(44100)
sink = PortAudioStream(0, 2; samplerate=sample_rate)

# write sine wave into sink only using PortAudio
# random frequencies, amplitude modulated
time = 0.1 # in seconds
reps = 50
freqs = 100 .+ 700*rand(reps) # in Hertz
for i in 1:reps
    x = cos.(2pi * (i-1) / reps) *
        sin.(2pi * (0:time*sample_rate-1) * freqs[i] / sample_rate)
    write(sink, x)
end
# frequency ramp
duration = 3
start, stop = 100, 300 # in Hertz
res = 100
restime = 0.1
rs = restime * sample_rate
freqs = range(start, stop, res)
x = [sin.(2pi * t * freqs[r]  / sample_rate)
    for r in 1:res for t in 0:rs]
write(sink, x)

# set up source and audio parameters
audio_channels = 2
source_rms = 0.5
source = NoiseSource(Float64, sample_rate, audio_channels, source_rms)

# run real time audio (static)
for frame = 1:100
    @pipe read(source, 0.01u"s") |> write(sink, _)
end

# run real time audio (modulated by amp)
amp = Amplification(current=0.1, target=0.3, change_limit=1e-3)
for frame = 1:250
    @pipe read(source, 0.01u"s") |> modify(amp, _) |> write(sink, _)
end

# run unbounded real time audio 
freq = 420
sine = SinusoidSource(Float64, 44100, freq)
is_open = true 
while is_open
    @pipe read(sine, 0.01u"s") |> modify(amp, _) |> write(sink, _)
end # interrupt to stop signal

# run frequency stack
freqs = 300:50:800
sine_stack = SinusoidSource(Float64, 44100, freqs)
for frame = 1:200
    @pipe read(sine_stack, 0.01u"s") |> modify(amp, _) |> write(sink, _)
end

# run frequency ramp
start = 100
stop = 420
num_frames = 200
for frame = 1:num_frames
    sine_stack = SinusoidSource(Float64, 44100,
        (1.0-frame/num_frames)*start + frame*stop/num_frames)
    @pipe read(sine_stack, 0.1u"s") |> modify(amp, _) |> write(sink, _)
end