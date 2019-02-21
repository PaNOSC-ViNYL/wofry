import numpy
from scipy.special import fresnel


def fraunhofer_analytical_rectangle(fresnel_number=None, propagation_distance=1.140, aperture_half=1e-3,
                                    wavelength=639e-9, detector_array=None, npoints=1000):

    if fresnel_number is None:
        fresnel_number = aperture_half**2 / (wavelength * propagation_distance)

    if detector_array is None:
        if fresnel_number > 1.0:
            window_aperture_ratio = 2.0
        else:
            window_aperture_ratio = 1.0 / fresnel_number

        x = numpy.linspace(-window_aperture_ratio * aperture_half,
                           window_aperture_ratio * aperture_half, npoints)
    else:
        x = detector_array.copy()

    argument_sinc = 2.0 * aperture_half * numpy.pi / wavelength / propagation_distance * x # TODO: check the 2??
    alpha = 2.0 * aperture_half / (wavelength * propagation_distance) ** (1.0 / 2.0) * \
            numpy.exp(1j * numpy.pi / wavelength / propagation_distance * x**2) * \
            numpy.sin(argument_sinc) / argument_sinc

    # TODO note that the global phase (Goldman 4-59) is missing

    return x, alpha


def fresnel_analytical_rectangle(fresnel_number=None, propagation_distance=1.140, aperture_half=1e-3, wavelength=639e-9,
                                 detector_array=None, npoints=1000):

    if fresnel_number is None:
        fresnel_number = aperture_half**2 / (wavelength * propagation_distance)

    if detector_array is None:
        if fresnel_number > 1.0:
            window_aperture_ratio = 2.0
        else:
            window_aperture_ratio = 1.0 / fresnel_number
        x = numpy.linspace(-window_aperture_ratio * aperture_half, window_aperture_ratio * aperture_half, npoints)
    else:
        x = detector_array.copy()

    s_plus  = numpy.sqrt(2.0 * fresnel_number) * ( 1.0 + x / aperture_half)
    s_minus = numpy.sqrt(2.0 * fresnel_number) * ( 1.0 - x / aperture_half)

    fs_plus,fc_plus = fresnel(s_plus)
    fs_minus,fc_minus = fresnel(s_minus)

    Ux = (fc_minus + fc_plus) + 1j * (fs_minus + fs_plus)
    Ux *= 1.0 / numpy.sqrt(2.0)

    # TODO note that the global phase (Goldman 4-59) is missing

    return x, Ux  # note that wavefield is being returned, not intensity!
