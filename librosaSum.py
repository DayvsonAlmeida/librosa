
y, sr = librosa.load('audio.wav', sr = 16000)

feat2 = librosa.feature.melspectrogram(y = y, sr = sr, n_mels = 20)

#melspectrogram
def melspectrogram(y=None, sr=22050, S=None, n_fft=2048, hop_length=512,
                   win_length=None, window='hann', center=True, pad_mode='reflect',
                   power=2.0, **kwargs):

    S = _spectrogram(y=y, n_fft=n_fft, hop_length=hop_length, power=power,
                            win_length=win_length, window=window, center=center,
                            pad_mode=pad_mode)

    # Build a Mel filter
    mel_basis = filters.mel(sr, n_fft, **kwargs)

    return np.dot(mel_basis, S) #produto

#_spectrogram
def _spectrogram(y=None, n_fft=2048, hop_length=512, power=2.0,
                 win_length=None, window='hann', center=True, pad_mode='reflect'):
        # Otherwise, compute a magnitude spectrogram from input
    return np.abs(stft(y, n_fft=n_fft, hop_length=hop_length,
                        win_length=win_length, center=center,
                        window=window, pad_mode=pad_mode))**power
    

def mel(sr=16000, n_fft=2048, n_mels=128, fmin=0.0, fmax=None, htk=False,
        norm=1, dtype=np.float32):

    if fmax is None:
        fmax = float(sr) / 2

    if norm is not None and norm != 1 and norm != np.inf:
        raise ParameterError('Unsupported norm: {}'.format(repr(norm)))

    # Initialize the weights
    n_mels = int(n_mels)
    weights = np.zeros((n_mels, int(1 + n_fft // 2)), dtype=dtype)

    # Center freqs of each FFT bin
    fftfreqs = fft_frequencies(sr=sr, n_fft=n_fft)                  #TO AQUI

    # 'Center freqs' of mel bands - uniformly spaced between limits
    mel_f = mel_frequencies(n_mels + 2, fmin=fmin, fmax=fmax, htk=htk)

    fdiff = np.diff(mel_f)
    ramps = np.subtract.outer(mel_f, fftfreqs)

    for i in range(n_mels):
        # lower and upper slopes for all bins
        lower = -ramps[i] / fdiff[i]
        upper = ramps[i+2] / fdiff[i+1]

        # .. then intersect them with each other and zero
        weights[i] = np.maximum(0, np.minimum(lower, upper))

    if norm == 1:
        # Slaney-style mel is scaled to be approx constant energy per channel
        enorm = 2.0 / (mel_f[2:n_mels+2] - mel_f[:n_mels])
        weights *= enorm[:, np.newaxis]

    # Only check weights if f_mel[0] is positive
    if not np.all((mel_f[:-2] == 0) | (weights.max(axis=1) > 0)):
        # This means we have an empty channel somewhere
        warnings.warn('Empty filters detected in mel frequency basis. '
                      'Some channels will produce empty responses. '
                      'Try increasing your sampling rate (and fmax) or '
                      'reducing n_mels.')

    return weights

#fft_frequencies
def fft_frequencies(sr=22050, n_fft=2048):
  
    return np.linspace(0,float(sr) / 2,int(1 + n_fft//2),endpoint=True)


def mel_frequencies(n_mels=128, fmin=0.0, fmax=11025.0, htk=False):

    min_mel = hz_to_mel(fmin, htk=htk)
    max_mel = hz_to_mel(fmax, htk=htk)

    mels = np.linspace(min_mel, max_mel, n_mels)

return mel_to_hz(mels, htk=htk)







