import numpy as np
import math
from scipy.interpolate import interp1d
import scipy.stats as stats
import matplotlib.pyplot as plt


def nanspectrum(TS, DT, TIME_UNITS, SEGMENTS, PLOT_OPTION, PLOT_BOOLEAN, INTERPMETHOD, *args):
    #Spectrum = np.fft.fft(TS) ** 2 / DT
    #f_vec = np.arange(1 / DT, 0.5 + 1 / DT, 1 / DT)
    #err = 0


    ## Alert the user if the time series cannot be properly segmented
    leng_TS = len(TS)
    #error_Msg = f"You need to make sure that your time series can actually be divided into {SEGMENTS} windows. As it is given, the time series has {leng_TS} data points. This can be factored into [{', '.join(map(str, factor(leng_TS)))}]."

    # Give yourself an error if you can't split TS into SEGMENTS segments:
    if leng_TS % SEGMENTS != 0:
        # Using one segment is not advised, because it will just result in the
        # absolute value of the fft, which will have a lot of noise, but if you
        # do need to use the full time series without segmenting, it would
        # probably be best to not dampen the edges and just use 'rectwin' for
        # the variable WINDOWMETHOD ("rectwin" gives all ones). A warning will be
        # given when WINDOWMETHOD is assigned.
        print(f"You need to make sure that your time series can actually be divided into the desired number of segments. Python doesn't seem to have a basic prime factoring function, which is very disappointing.")
        raise ValueError("The time series is of the wrong length.")
    # else:
        # Display the warning message
        # print(error_Msg)

    ## Handle WINDOWMETHOD variable:
    N = leng_TS / SEGMENTS

    if len(args) == 0:
        if SEGMENTS == 1:
            WINDOWMETHOD = 'np.hanning'
            # should be "rectwin" but numpy doesn't have that
            # possible options include: bartlett,blackman,hamming,hanning,kaiser
        else:
            WINDOWMETHOD = 'np.hanning'
    elif len(args) == 1:
        WINDOWMETHOD = args[0]
        if isinstance(WINDOWMETHOD, str):
            if SEGMENTS == 1:
                if WINDOWMETHOD not in ['rectwin', 'boxcar']:
                    print('You have chosen not to segment your time series. In this case, it is advised (but not required) that you choose "rectwin" or "boxcar" for WINDOWMETHOD, because they do not deemphasize the edges of the time series.')
            else:
                pass
        elif isinstance(WINDOWMETHOD, list) and len(WINDOWMETHOD) == N:
            WINDOWMETHOD = WINDOWMETHOD
        elif isinstance(WINDOWMETHOD, list) and len(WINDOWMETHOD) == N:
            WINDOWMETHOD = WINDOWMETHOD
        else:
            raise ValueError('The optional eighth argument "WINDOWMETHOD" needs to be a string (starting with "np.") for a window function or a vector of length len(TS)/SEGMENTS.')
    else:
        raise ValueError('"nanspectrum" only takes 7 or 8 arguments.')

    ## Make a column into a row
    if TS.ndim == 1 and TS.shape[0] > 1:  # Check if TS is a row vector
        pass  # No action needed if TS is a row vector
    else:  # Reformat if given as a column
        TS = np.transpose(TS)

    ## Interpolate over NaN's using the chosen method "INTERPMETHOD"
    T = np.arange(1, leng_TS + 1)
    if isinstance(INTERPMETHOD, str):
        finite_indices = np.isfinite(TS)
        interp_indices = np.isfinite(T)
        interp_func = interp1d(T[finite_indices], TS[finite_indices], kind=INTERPMETHOD)
        TS = interp_func(T[interp_indices])
    elif INTERPMETHOD == 0:
        TS_nonan = TS[np.isfinite(TS)]
        #$ print(T.shape) #$ print(TS.shape)
        T_nonan = T[np.isfinite(TS)]
        A = np.ones((len(TS_nonan), 2))
        A[:, 1] = T_nonan
        m = np.linalg.lstsq(A, TS_nonan, rcond=None)[0]
        nancoords = np.where(~np.isfinite(TS))[0]
        for i in nancoords:
            TS[i] = m[0] + m[1] * T[i]
    elif INTERPMETHOD == 1:
        TS_nonan = TS[np.isfinite(TS)]
        T_nonan = T[np.isfinite(TS)]
        A = np.ones((len(TS_nonan), 2))
        A[:, 1] = T_nonan
        m = np.linalg.lstsq(A, TS_nonan, rcond=None)[0]
        TS_nonandetrended = TS_nonan - np.dot(A, m)
        TS_nandetrended = TS.copy()
        TS_nandetrended[np.isfinite(TS)] = TS_nonandetrended
        TS_var = np.nanvar(TS_nandetrended)
        nancoords = np.where(~np.isfinite(TS))[0]
        for i in nancoords:
            TS[i] = np.sqrt(TS_var) * np.random.randn() + m[0] + m[1] * T[i]
    else:
        help(nanspectrum)
        raise ValueError('The variable "INTERPMETHOD" must be 0, 1, or one of several specific strings. See the documentation above.')

    ## Segment (reduce noise), window (account for edge effects), and calculated "Spectrum"
    N = len(TS) // SEGMENTS

    if isinstance(WINDOWMETHOD, str) and (WINDOWMETHOD not in ['rectwin','boxcar']):
        Window = eval("np.tile(" + WINDOWMETHOD + "(N), (1, 2*SEGMENTS - 1))")
    elif isinstance(WINDOWMETHOD, str) and (WINDOWMETHOD in ['rectwin','boxcar']):
        Window = np.tile(np.hanning(N), (1, 2*SEGMENTS - 1))
        Window = Window*0 + 1 # <--- equivalent to rectwin and boxcar
    else:
        Window = np.tile(WINDOWMETHOD, (1, 2*SEGMENTS - 1))
    Window = Window[0]
    #plt.plot(Window)#$

    # norm_factor = 1 / np.mean([0, Window[:, 0]]**2)
    norm_factor = 1 / np.mean(([0] + Window[:])**2)

    seg_TS = np.zeros((N, 2*SEGMENTS - 1))
    for i in range(0, 2*SEGMENTS - 1, 2):
        seg_TS[:, i] = TS[((i//2)*N):(((i+2)//2)*N)]

    Offset = N // 2
    for i in range(1, 2*SEGMENTS - 1, 2):
        #seg_TS[:, i] = TS[(Offset + 1 + ((i-1)//2)*N):(Offset + ((i+1)//2)*N)]
        seg_TS[:, i] = TS[(Offset + (i//2)*N):(Offset + ((i+2)//2)*N)]

    seg_TS = seg_TS - np.mean(seg_TS, axis=0)
    window_seg_TS = np.reshape(Window, seg_TS.shape, order='F') * seg_TS

    fft_window_seg_TS = np.fft.fft(window_seg_TS, axis=0)

    if N % 2 == 0:
        amp = np.abs(fft_window_seg_TS[0:(N//2 + 1), :])**2
    elif N % 2 == 1:
        amp = np.abs(fft_window_seg_TS[0:((N+1)//2), :])**2
    else:
        raise ValueError('Fatal error, figure it out yourself')

    if N % 2 == 0:
        amp[1:(-1), :] = 2 * amp[1:(-1), :]
    elif N % 2 == 1:
        amp[1:, :] = 2 * amp[1:, :]
    else:
        raise ValueError('Fatal error, figure it out yourself')

    f_s = 1 / DT
    timelength_record = ((len(TS) * DT) / SEGMENTS)

    f_Ny = 0.5 * f_s
    df = 1 / timelength_record

    if N % 2 == 0:
        f_vec = np.arange(df, f_Ny + df, df)
    elif N % 2 == 1:
        f_vec = np.arange(df, f_Ny - df/2 + df, df)
    else:
        raise ValueError('Fatal error, figure it out yourself')

    amp = norm_factor * amp * DT / N
    Spectrum = np.mean(amp[1:, :], axis=1)

    ## Error
    # In short, do this calculation for 95% confidence ratio:
    EffectiveSegments = 2 * SEGMENTS - 1
    # Degrees of freedom; for overlapping Hann windowed segments,
    # this is the total number of
    # segments, i.e.:
    # 2*SEGMENTS - 1
    # This may only work perfectly with the Hann window.
    # With other windows, it's not entirely clear, but my
    # guess is, if the window is "Hann-like", it will not
    # make a big difference.

    err_high = 2 * EffectiveSegments / stats.chi2.ppf(0.05 / 2, 2 * EffectiveSegments)
    err_low = 2 * EffectiveSegments / stats.chi2.ppf(1 - 0.05 / 2, 2 * EffectiveSegments)
    err = [err_low, err_high]

    ## Plot
    if PLOT_BOOLEAN:
        # Sort out blank plot option
        if PLOT_OPTION == 0:
            PLOT_OPTION = '.-'
        else:
            pass

        # One can comment out "plt.figure()" to facilitate plotting
        # several spectra together from calling this function multiple separate times.
        plt.figure()
        plt.loglog(f_vec, Spectrum, PLOT_OPTION)
        #plt.plot(f_vec, Spectrum, PLOT_OPTION)
        plt.xlabel('Cycles per ' + TIME_UNITS)#, interpreter='latex'
        plt.ylabel('Spectral density [(time\_series\_units)$^2$ (cycles per ' +
                   TIME_UNITS + ')$^{-1}$]')#, interpreter='latex'

        # Error bar plotted
        plt.loglog(np.array([1.1, 1.1]) * f_vec[-1], np.array([1, err_high / err_low]) * min(Spectrum), 'k')

        #plt.ylim((80,10**3))
        plt.show()
    else:
        pass

    return Spectrum, f_vec, err
