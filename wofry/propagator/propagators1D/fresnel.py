import numpy

from srxraylib.util.data_structures import ScaledArray
from wofry.propagator.generic_wavefront import GenericWavefront1D
from wofry.propagator.propagator import Generic1DPropagator

class Fresnel1D(Generic1DPropagator):

    def get_handler_name(self):
        return "FRESNEL_1D"

    """
    1D Fresnel propagator using convolution via Fourier transform
    :param wavefront:
    :param propagation_distance: propagation distance
    :return: a new 1D wavefront object with propagated wavefront
    """
    def do_specific_progation(self, wavefront, propagation_distance, parameters):

        fft_scale = numpy.fft.fftfreq(wavefront.size())/wavefront.delta()

        fft = numpy.fft.fft(wavefront.get_complex_amplitude())
        fft *= numpy.exp((-1.0j) * numpy.pi * wavefront.get_wavelength() * propagation_distance * fft_scale**2)
        ifft = numpy.fft.ifft(fft)

        return GenericWavefront1D(wavefront.get_wavelength(), ScaledArray.initialize_from_steps(ifft, wavefront.offset(), wavefront.delta()))

class FresnelConvolution1D(Generic1DPropagator):

    def get_handler_name(self):
        return "FRESNEL_CONVOLUTION_1D"

    """
    1D Fresnel propagator using direct convolution
    :param wavefront:
    :param propagation_distance:
    :return:
    """
    def do_specific_progation(self, wavefront, propagation_distance, parameters):
        # instead of numpy.convolve, this can be used:
        # from scipy.signal import fftconvolve

        kernel = numpy.exp(1j*2*numpy.pi/wavefront.get_wavelength() * wavefront.get_abscissas()**2 / 2 / propagation_distance)
        kernel *= numpy.exp(1j*2*numpy.pi/wavefront.get_wavelength() * propagation_distance)
        kernel /=  1j * wavefront.get_wavelength() * propagation_distance
        tmp = numpy.convolve(wavefront.get_complex_amplitude(),kernel,mode='same')

        return GenericWavefront1D(wavefront.get_wavelength(), ScaledArray.initialize_from_steps(tmp, wavefront.offset(), wavefront.delta()))