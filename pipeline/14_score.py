"""Original ambient score for the combined video, synthesised from scratch with numpy (no samples).

Structure follows the video's scenes (seconds are passed in by the caller, defaults below):
  pads    slow chord progression D maj9 - B m11 - G maj7(#11) - A sus2, soft detuned saw-ish voices
  bass    chord roots, sine with gentle saturation
  swell   filtered noise + rising pad during the flight to the star
  sparkle pentatonic plinks over the starfield sections
  pulse   one soft low "heartbeat" each time the planet crosses the star
A cheap stereo reverb (parallel feedback combs + allpass) glues it together.

Writes results/media/score.wav (44.1 kHz stereo). Usage: python pipeline/14_score.py [DURATION_S]
"""
import sys
import wave

import env  # noqa: F401

import numpy as np
from scipy.signal import lfilter

from env import RESULTS

SR = 44100
RNG = np.random.default_rng(370612437)
A4 = 440.0
# chord voicings as semitones from A4 (D maj9, B m11, G maj7#11, A sus2), each held CHORD_S seconds
CHORDS = [[-19, -12, -8, -5, -1, 2], [-22, -15, -10, -7, -3, 0], [-26, -19, -14, -10, -7, -3], [-24, -12, -10, -5, 0, 2]]
ROOTS = [-31, -34, -38, -36]
CHORD_S = 8.0
PENTATONIC = [0, 2, 4, 7, 9]  # D major pentatonic over the D tonic (semitones above D5)
D5 = -7  # semitones from A4


def hz(semi):
    return A4 * 2 ** (semi / 12)


def env_adsr(n, attack, release):
    e = np.ones(n)
    a, r = min(int(attack * SR), n // 2), min(int(release * SR), n // 2)
    e[:a] = np.linspace(0, 1, a) ** 2
    e[-r:] *= np.linspace(1, 0, r) ** 2
    return e


def lowpass(x, cutoff, block=512):
    """One-pole low-pass; cutoff may vary with time (held constant within 512-sample blocks)."""
    cutoff = np.broadcast_to(cutoff, x.shape)
    y = np.empty_like(x)
    zi = np.zeros(1)
    for i in range(0, len(x), block):
        alpha = 1 - np.exp(-2 * np.pi * cutoff[i] / SR)
        y[i:i + block], zi = lfilter([alpha], [1, alpha - 1], x[i:i + block], zi=zi)
    return y


def pad_voice(f, n):
    t = np.arange(n) / SR
    out = np.zeros(n)
    for detune in (-0.07, 0.0, 0.06):  # three slightly detuned oscillators -> chorus shimmer
        ff = f * 2 ** (detune / 12)
        ph = RNG.uniform(0, 2 * np.pi)
        for k in range(1, 7):  # band-limited saw-ish: harmonics falling as 1/k^1.6
            out += np.sin(2 * np.pi * ff * k * t + ph * k) / k ** 1.6
    return out / 3


def reverb(x, mix=0.3):
    """Parallel feedback combs (different delays per ear) for a wide, soft room."""
    out = np.zeros((len(x), 2))
    for ch, delays in enumerate(([1557, 1617, 1491, 1422], [1277, 1356, 1188, 1116])):
        wet = np.zeros(len(x))
        for d in delays:
            y = x.copy()
            for i in range(d, len(x), d):  # block-wise recursion y[n] += g * y[n - d]
                y[i:i + d] += 0.8 * y[i - d:i][: len(y[i:i + d])]
            wet += y
        wet = lowpass(wet / len(delays), 3500)
        out[:, ch] = (1 - mix) * x + mix * wet * (np.std(x) / (np.std(wet) + 1e-12))
    return out


def score(duration, flight=(10.0, 26.0), transits=(), sparkle_spans=((4.0, 26.0), (46.0, 62.0))):
    n = int(duration * SR)
    t = np.arange(n) / SR
    mono = np.zeros(n)
    # pads and bass, chords overlapping by 2 s
    start = 0.0
    i = 0
    while start < duration:
        seg = int(min(CHORD_S + 3, duration - start + 0.01) * SR)
        s0 = int(start * SR)
        seg = min(seg, n - s0)
        if seg <= 0:
            break
        chord = CHORDS[i % 4]
        e = env_adsr(seg, 2.5, 3.0)
        mono[s0:s0 + seg] += 0.09 * e * sum(pad_voice(hz(s), seg) for s in chord)
        bass = np.sin(2 * np.pi * hz(ROOTS[i % 4]) * np.arange(seg) / SR)
        mono[s0:s0 + seg] += 0.22 * e * np.tanh(1.5 * bass)
        start += CHORD_S
        i += 1
    # darker filter early, opening up during the flight
    f0, f1 = flight
    cutoff = 900 + 3200 * np.clip((t - f0) / (f1 - f0), 0, 1) ** 1.5
    cutoff[t > f1] = 4100 - 1800 * np.clip((t[t > f1] - f1) / 6, 0, 1)
    mono = lowpass(mono, cutoff)
    # flight swell: band-limited noise rising, then cut at arrival
    sw = (t > f0) & (t < f1 + 1.5)
    noise = lowpass(lowpass(RNG.standard_normal(n), 200 + 1300 * np.clip((t - f0) / (f1 - f0), 0, 1)), 1800)
    swell = np.clip((t - f0) / (f1 - f0), 0, 1) ** 2.2 * np.exp(-np.clip(t - f1, 0, None) * 1.5)
    mono += 0.12 * noise * swell * sw
    # sparkle plinks
    plinks = np.zeros(n)
    for a, b in sparkle_spans:
        for tp in np.sort(RNG.uniform(a, b, int((b - a) * 1.3))):
            semi = D5 + 12 + RNG.choice(PENTATONIC) + 12 * RNG.integers(0, 2)
            k0 = int(tp * SR)
            m = min(int(2.5 * SR), n - k0)
            tt = np.arange(m) / SR
            plinks[k0:k0 + m] += 0.05 * np.sin(2 * np.pi * hz(semi) * tt) * np.exp(-tt * 2.2) * (1 - np.exp(-tt * 300))
    mono += plinks
    # soft low pulse at each transit mid-time
    for tp in transits:
        k0 = int(tp * SR)
        m = min(int(1.8 * SR), n - k0)
        tt = np.arange(m) / SR
        mono[k0:k0 + m] += 0.35 * np.sin(2 * np.pi * 55 * tt) * np.exp(-tt * 2.5) * (1 - np.exp(-tt * 60))
    st = reverb(mono)
    # master fades (stereo width comes from the reverb; a Haas delay would cancel the bass in mono)
    fade = np.clip(t / 2.0, 0, 1) * np.clip((duration - t) / 4.0, 0, 1)
    st *= fade[:, None]
    return st / np.max(np.abs(st)) * 0.7


def write_wav(path, stereo):
    data = (np.clip(stereo, -1, 1) * 32767).astype("<i2")
    with wave.open(str(path), "wb") as w:
        w.setnchannels(2)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes(data.tobytes())


def main():
    duration = float(sys.argv[1]) if len(sys.argv) > 1 else 66.0
    out = RESULTS / "media"
    out.mkdir(parents=True, exist_ok=True)
    write_wav(out / "score.wav", score(duration))
    print(f"wrote {out / 'score.wav'} ({duration:.0f} s)")


if __name__ == "__main__":
    main()
