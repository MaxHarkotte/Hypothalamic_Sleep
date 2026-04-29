import jax
import jax.numpy as jnp
import xarray as xr
import numpy as np


def mt_specpb_xr_chunked_jax(
    data, Fs=1000, NW=4, trial_length=None, has_trial_dim=False
):
    """
    Multitaper spectral estimate supporting:
        - data(channel, time) with chunking
        - data(channel, trial, time) without chunking

    Output dims: (channel, trial, freq)
    """

    # -------------------------
    # Case 1: (C,Tr,T) already has trials
    # -------------------------
    if has_trial_dim:
        if (
            ("channel" not in data.dims)
            or ("trial" not in data.dims)
            or ("time" not in data.dims)
        ):
            raise ValueError(
                "Data must have dims ('channel','trial','time') when has_trial_dim=True."
            )

        x = jnp.asarray(data.values)  # (C,Tr,Ttrial)
        nchan, ntrials, trial_length = x.shape

    # -------------------------
    # Case 2: (C,T) → chunk trials
    # -------------------------
    else:
        if ("channel" not in data.dims) or ("time" not in data.dims):
            raise ValueError(
                "Data must have dims ('channel','time') when has_trial_dim=False."
            )

        x = jnp.asarray(data.values)  # (C,T)
        nchan, ntime = x.shape

        # Default: one trial spanning entire recording
        if trial_length is None:
            trial_length = ntime

        ntrials = ntime // trial_length
        if ntrials == 0:
            raise ValueError("trial_length > signal length")

        # Truncate and reshape into trials: (C,Tr,Ttrial)
        x = x[:, : ntrials * trial_length]
        x = x.reshape(nchan, ntrials, trial_length)

    tapers, _ = dpss_windows(trial_length, NW, 2 * NW - 1)
    tapers = jnp.asarray(tapers) * jnp.sqrt(Fs)  # (K,Ttrial)
    k = tapers.shape[0]
    T = jnp.fft.rfft(tapers, axis=-1)
    freqs = jnp.fft.rfftfreq(trial_length, 1 / Fs)

    # -------------------------
    # JAX FFT computation
    # -------------------------
    @jax.jit
    def compute(x):
        """
        x: (C,Tr,Ttrial)
        output: (C,Tr,nf)
        """

        xt = (
            x[:, :, None, :] * tapers[None, None, :, :]
        )  # Broadcast multiply → (C,Tr,K,Ttrial)

        J = jnp.fft.rfft(xt, axis=-1)  # (C,Tr,K,nf)

        # DC term per channel/per trial
        mean_val = jnp.mean(x, axis=-1)  # (C,Tr)
        dc = mean_val[:, :, None, None] * T[None, None, :, :]  # (C,Tr,K,nf)

        P = jnp.real((J - dc) * jnp.conj(J - dc))  # power

        return jnp.mean(P, axis=2)  # average over tapers → (C,Tr,nf)

    S = compute(x)

    # -------------------------
    # Return DataArray
    # -------------------------
    return xr.DataArray(
        S,
        dims=("channel", "trial", "freq"),
        coords={
            "channel": data.channel,
            "trial": jnp.arange(ntrials),
            "freq": freqs,
        },
        attrs=data.attrs,
    )


def mt_specpb(data, Fs=1000, NW=4, chunk_size=None, chunk_avg=False):
    tapers, _ = dpss_windows(data.shape[-1], NW, 2 * NW - 1)  # Compute the tapers,
    tapers *= np.sqrt(Fs)  # ... and scale them.
    if chunk_size is None:
        chunk_size = data.shape[0]

    nchunks = int(np.ceil(data.shape[0] / chunk_size))
    spectra = []
    spectra_sem = []

    for i in range(nchunks):
        chunk = data[i * chunk_size : (i + 1) * chunk_size, :]
        # Taper and FFT
        dataT = np.array(
            [[trial * t for t in tapers] for trial in chunk]
        )  # shape: (nchunk, k, time)
        T = np.fft.rfft(tapers, axis=-1)  # shape: (k, nf)
        J = np.fft.rfft(dataT, axis=-1)  # shape: (nchunk, k, nf)

        # Subtract DC
        dc = np.array([T * trial.mean() for trial in chunk])  # shape: (nchunk, k, nf)
        J -= dc

        # Spectrum: power
        J *= J.conj()  # power
        S_chunk = J.mean(1).real
        spectra.append(np.mean(S_chunk, axis=0))  # mean across chunk
        spectra_sem.append(stats.sem(S_chunk, axis=0))
    spectra = np.stack(spectra)  # shape: (nchunks, nf)
    spectra_sem = np.stack(spectra_sem)
    f = np.fft.rfftfreq(data.shape[-1], 1 / Fs)
    if chunk_avg:
        spectra = spectra.mean(0)  # Average across trials.
        spectra_sem = spectra_sem.mean(0)
    return f, spectra, spectra_sem


def mt_specpb_xr_chunked_jax_dep(data, Fs=1000, NW=4, trial_length=None):
    """
    Multitaper spectral estimate with trial chunking.
    Input dims:  (channel, time)
    Output dims: (channel, trial, freq)
    """

    if ("channel" not in data.dims) or ("time" not in data.dims):
        raise ValueError("Data must have dims ('channel','time').")
    x = jnp.asarray(data.values)  # (C,T)
    nchan, ntime = x.shape
    if trial_length is None:
        trial_length = ntime

    ntrials = ntime // trial_length
    if ntrials == 0:
        raise ValueError("trial_length > signal length")

    # Truncate and reshape into trials: (C, ntrials, Ttrial)
    x = x[:, : ntrials * trial_length]
    x = x.reshape(nchan, ntrials, trial_length)
    tapers, _ = dpss_windows(trial_length, NW, 2 * NW - 1)
    tapers = jnp.asarray(tapers) * jnp.sqrt(Fs)  # (K,Ttrial)
    k = tapers.shape[0]
    T = jnp.fft.rfft(tapers, axis=-1)
    freqs = jnp.fft.rfftfreq(trial_length, 1 / Fs)

    # ---- JAX computation --------------------------------------
    @jax.jit
    def compute(x):
        """
        x: (C,Tr,Ttrial)
        output: (C,Tr,nf)
        """

        xt = (
            x[:, :, None, :] * tapers[None, None, :, :]
        )  # Broadcast multiply: (C,Tr,1,Ttrial) * (1,1,K,Ttrial)
        J = jnp.fft.rfft(xt, axis=-1)  # FFT: (C,Tr,K,nf)

        mean_val = jnp.mean(x, axis=-1)  # (C,Tr) # DC term per channel/per trial
        dc = mean_val[:, :, None, None] * T[None, None, :, :]  # (C,Tr,K,nf)
        P = jnp.real(
            (J - dc) * jnp.conj(J - dc)
        )  # (C,Tr,K,nf) # Subtract DC and compute power
        return jnp.mean(P, axis=2)  # → (C,Tr,nf) # Average over tapers

    S = compute(x)
    return xr.DataArray(
        S,
        dims=("channel", "trial", "freq"),
        coords={
            "channel": data.channel,
            "trial": jnp.arange(ntrials),
            "freq": freqs,
        },
        attrs=data.attrs,
    )
