# -*- coding: utf-8 -*-
#from __future__ import print_function

"""
Created on Mon Apr  4 15:32:35 2016

@author: fwinkler
"""

import numpy as np
import re
import os


class WavimgPrm(object):
    """
    Wavimg object that can be used to load, modify and save a wavimg parameter
    file for the Dr. Probe software package.

    Attributes
    ----------
    wave_files : string
        Wave function file name string used to locate existing wave functions. Use quotation
        marks when the string includes space characters.
        Default : 'wav/xxx.wav'
    wave_dim : (int, int)
        Dimension of the wave data in pixels, <nx> = number of horizontal wave pixels,
        <ny> = number of vertical wave pixels.
        Default : (0, 0)
    wave_sampling : (float, float)
        Sampling rate of the wave data (<sx> = horizontal, <sy> = vertical) [nm/pix].
        Default : (0, 0)
    high_tension : int
        TEM high-tension used for wave function calculation [kV].
        Default : 80
    output_format : int
        Image output type option:
        0 = TEM image,
        1 = complex image plane wave,
        2 = wave amplitude,
        3 = wave phase,
        4 = wave real part,
        5 = wave imaginary part,
        6 = TEM image map of 2 variables
        Default : 1
    output_files : string
        Image output file name string. Use quotation marks when the string includes space
        characters.
        Default : 'img/xxx.dat'
    output_dim : (int, int)
        Image output size (<ix> = horizontal , <iy> = vertical) in number of pixels.
        Default : (0, 0)
    noise : (int, int, int, int)
        Flag and parameters for creating integer images with optional noise. Flag <intflg> 0 =
        off (default), 1 = 32-bit, 2 = 16-bit,
        Parameter:
        <mean> = mean vacuum intensity,
        <conv> = electron to counts conversion rate,
        <rnoise> detector readout rms noise level in counts.
        Default : (0, 1, 1, 0)
    flag_spec_frame : int
        Flag activating the extraction of a special image frame (0=OFF, 1=ON). The frame
        parameters are defined in the lines below.
        Default : 0
    output_sampling : float
        Image output sampling rate [nm/pix], isotropic. The parameter is used only if
        'flag_spec_frame' is set to 1.
        Default : 0
    img_frame_offset : (int, int)
        Image frame offset in pixels of the input wave. The parameter is used only if
        'flag_spec_frame' is set to 1.
        Default : (0, 0)
    img_rot : int
        Image frame rotation in [deg] with respect to the input wave horizontal axis. The
        parameter is used only if 'flag_spec_frame' is set to 1.
        Default : 0
    coherence_model : int
        Coherence calculation model switch:
        1 = averaging of coherent sub images explicit focal variation but
            quasi-coherent spatial envelope,
        2 = averaging of coherent sub images with explicit focal and angular
            variation,
        3 = quasi-coherent linear envelopes,
        4 = Fourier-space synthesis with partially coherent TCC,
        5 =  averaging of coherent sub images with explicit focal, angular,
            and frozen lattice variation)
        Default : 1
    temp_coherence : (int, float)
        Flag and parameters for partial temporal coherence: <ptcflg> = flag (0=OFF, 1=ON),
        <f-spread> = focus spread (1/e) half width [nm].
        Default : (1, 0.5)
    spat_coherence : (int, float)
        Flag and parameters for partial spatial coherence: <pscflg> = flag (0=OFF, 1=ON),
        <s-conv> = beam convergence (1/e) half width [mrad].
        Default : (1, 0.4)
    mtf : (int, float, string)
        Flag and parameters for applying the detector MTF: <mtfflag> = flag (0=OFF, 1=ON),
        <mtf-scale> = calculation scale of the mtf = (sampling rate experiment)/(sampling rate
        simulation), <mtf-file> = File name string to locate the MTF data. Use quotation marks
        when the string includes space characters.
        Default : (1, 1, 'PICO-US4k-080_mtf_bin1_4096.mtf')
    vibration : (int, float, float, int)
        Flag and parameters for a vibration envelope:
        <vibflg> = flag (0=OFF, 1=ON-ISO, 2=ON-ANISO),
        <vibprm1>, <vibprm1> = vibration RMS amplitudes [nm],
        <vibprm3> = orientation [deg] of the primary vibration amplitude w.r.t. the horizontal
        image axis.
        Default : (1, 0.022, 0.022, 0)
    number_of_aber : int
        Number of aberrations defined
        Default : 0
    oa_radius : int
        Objective aperture radius [mrad]. Set to very large values to deactivate.
        Default = 15
    oa_position : (int, int)
        Center of the objective aperture with respect to the zero beam [mrad].
        Default = (0, 0)

    Parameters
    ----------
    wavimg_dict : dictionary
        A dictionary with all attributes written to the parameter file. Can be loaded from an
        existing parameter file.
    aberrations_dict : dictionary
        A dictionary with all aberrations that will be included. Only those aberrations,
        whose values are non zero count. Can be loaded from an existing parameter file.
        Aberration indices:
            0: 'image_shift',
            1: 'defocus',
            2: '2-fold-astigmatism',
            3: 'coma',
            4: '3-fold-astigmatism',
            5: 'CS',
            6: 'star_aberration',
            7: '4-fold-astigmatism',
            8: 'coma(5th)',
            9: 'lobe-aberration',
            10: '5-fold-astigmatism',
            11: 'C5'

    Methods
    -------
    load_wavimg_prm(prm_filename)
        Creates a wavimgprm object with all attributes taken from the parameterfile
        'prm_filename'.
    save_wavimg_prm(prm_filename)
        Saves a wavimgprm object to a parameterfile. The path and name of the parameterfile is
        defined by prm_filename.

    """

    def __init__(self, wavimg_dict=None, aberrations_dict=None):
        if wavimg_dict is None:
            wavimg_dict = {}

        if aberrations_dict is None:
            self.aberrations_dict = {}
        else:
            self.aberrations_dict = aberrations_dict

        self.wave_files = wavimg_dict.get('wave_files', "wav/xxx.wav")
        self.wave_dim = wavimg_dict.get('wave_dim', (0, 0))
        self.wave_sampling = wavimg_dict.get('wave_sampling', (0, 0))
        self.high_tension = wavimg_dict.get('high_tension', 80)
        self.output_format = wavimg_dict.get('output_format', 1)
        self.output_files = wavimg_dict.get('output_files', "img/xxx.dat")
        self.output_dim = wavimg_dict.get('output_dim', (0, 0))
        self.noise = wavimg_dict.get('noise', (0, 1, 1, 0))
        self.flag_spec_frame = wavimg_dict.get('flag_spec_frame', 0)
        self.output_sampling = wavimg_dict.get('output_sampling', 0)
        self.img_frame_offset = wavimg_dict.get('img_frame_offset', (0, 0))
        self.img_rot = wavimg_dict.get('img_rot', 0)
        self.coherence_model = wavimg_dict.get('coherence_model', 1)
        self.temp_coherence = wavimg_dict.get('temp_coherence', (1, 0.5))
        self.spat_coherence = wavimg_dict.get('spat_coherence', (1, 0.4))
        self.mtf = wavimg_dict.get('mtf', (1, 1, "PICO-US4k-080_mtf_bin1_4096.mtf"))
        self.vibration = wavimg_dict.get('vibration', (1, 0.022, 0.022, 0))
        # self.number_of_aber = wavimg_dict.get('number_of_aber', self.number_of_aberrations)
        self.oa_radius = wavimg_dict.get('oa_radius', 15)
        self.oa_position = wavimg_dict.get('oa_position', (0, 0))
        self.number_of_loops = wavimg_dict.get('number_of_loops', 0)

    @property
    def number_of_aberrations(self):
        return len(self.aberrations_dict.keys())

    def load_wavimg_prm(self, prm_filename, output=False):
        """
        Loads the parameterfile 'prm_filename' .

        Parameters
        ----------
        prm_filename : str
            The name of the parameterfile.
        output : bool, optional
            Flag for terminal output.
        """
        aberrations_dict = {}

        with open(prm_filename, 'r') as prm:
            # rad content from prm file and split at comma and/or spaces
            content = prm.readlines()
            content = [re.split(r'[,\s]\s*', line) for line in content]
            # read lines to dictionary wavimg_dict
            self.wave_files = content[0][0]
            self.wave_dim = (int(content[1][0]), int(content[1][1]))
            self.wave_sampling = (float(content[2][0]), float(content[2][1]))
            self.high_tension = int(content[3][0])
            self.output_format = int(content[4][0])
            self.output_files = content[5][0]
            self.output_dim = (int(content[6][0]), int(content[6][1]))
            self.noise = (int(content[7][0]), int(content[7][1]), int(content[7][2]),
                          int(content[7][3]))
            self.flag_spec_frame = int(content[8][0])
            self.output_sampling = float(content[9][0])
            self.img_frame_offset = (int(content[10][0]), int(content[10][1]))
            self.img_rot = float(content[11][0])
            self.coherence_model = int(content[12][0])
            self.temp_coherence = (int(content[13][0]), float(content[13][1]))
            self.spat_coherence = (int(content[14][0]), float(content[14][1]))
            self.mtf = (int(content[15][0]), float(content[15][1]), content[15][2])
            self.vibration = (int(content[16][0]), float(content[16][1]), float(content[16][2]),
                              float(content[16][3]))
            n = int(content[17][0])
            # self.number_of_aber = n
            for i in range(n):
                index = int(content[18 + i][0])
                aberrations_dict[index] = (float(content[18 + i][1]), float(content[18 + i][2]))
            #self.number_of_aber = int(content[17][0])
            self.aberrations_dict = aberrations_dict
            self.oa_radius = float(content[18 + n][0])
            self.oa_position = (int(content[19 + n][0]), int(content[19 + n][1]))
            self.number_of_loops = int(content[20 + n][0])

        if output:
            print("Parameters successfully loaded from file '{}'!".format(prm_filename))

    def save_wavimg_prm(self, prm_filename, output=False):
        """
        Saves the WavimgPrm object in the parameterfile 'prm_filename'.

        Parameters
        ----------
        prm_filename : str
            The path and name of the parameterfile.
        output : bool, optional
            Flag for terminal output
        """

        # directory = prm_filename.rsplit('/', 1)[0]
        directory = os.path.split(prm_filename)[0]
        if directory:
            if not os.path.isdir(directory):
                os.makedirs(directory, exist_ok=True)

        aberrations = {0: 'image_shift',
                       1: 'defocus',
                       2: '2-fold-astigmatism',
                       3: 'coma',
                       4: '3-fold-astigmatism',
                       5: 'CS',
                       6: 'star_aberration',
                       7: '4-fold-astigmatism',
                       8: 'coma(5th)',
                       9: 'lobe-aberration',
                       10: '5-fold-astigmatism',
                       11: 'C5'}

        with open(prm_filename, 'w') as prm:
            string_1 = "Wave function file name string used to locate existing wave functions. " \
                       "Use quotation marks when the string includes space characters."
            prm.write("'{}' ! {}\n".format(self.wave_files, string_1))
            string_2 = "Dimension of the wave data in pixels, <nx> = number of horizontal wave " \
                       "pixels, <ny> = number of vertical wave pixels."
            prm.write("{}, {} ! {}\n".format(self.wave_dim[0], self.wave_dim[1], string_2))
            string_3 = "Sampling rate of the wave data (<sx> = horizontal, <sy> = vertical) [" \
                       "nm/pix]."
            prm.write("{}, {} ! {}\n".format(self.wave_sampling[0], self.wave_sampling[1],
                                             string_3))
            string_4 = "TEM high-tension used for wave function calculation [kV]."
            prm.write("{} ! {}\n".format(self.high_tension, string_4))
            string_5 = "Image output type option: 0 = TEM image, 1 = complex image plane wave, " \
                       "2 = wave amplitude, 3 = wave phase, 4 = wave real part, 5 = wave " \
                       "imaginary part, 6 = TEM image map of 2 variables."
            prm.write("{} ! {}\n".format(self.output_format, string_5))
            string_6 = "Image output file name string. Use quotation marks when the string " \
                       "includes space characters."
            prm.write("'{}' ! {}\n".format(self.output_files, string_6))
            string_7 = "Image output size (<ix> = horizontal , <iy> = vertical) in number of " \
                       "pixels."
            prm.write("{}, {} ! {}\n".format(self.output_dim[0], self.output_dim[1], string_7))
            string_8 = "Flag and parameters for creating integer images with optional noise. " \
                       "Flag <intflg> 0 = off (default), 1 = 32-bit, 2 = 16-bit, Parameter: " \
                       "<mean> = mean vacuum intensity, <conv> = electron to counts conversion " \
                       "rate, <rnoise> detector readout rms noise level in counts."
            prm.write("{}, {}, {}, {} ! {}\n".format(self.noise[0], self.noise[1], self.noise[2],
                                                     self.noise[3], string_8))
            string_9 = "Flag activating the extraction of a special image frame (0=OFF, " \
                       "1=ON). The frame parameters are defined in the lines below."
            prm.write("{} ! {}\n".format(self.flag_spec_frame, string_9))
            string_10 = "Image output sampling rate [nm/pix], isotropic. The parameter is used " \
                        "only if the Flag in line 09 is set to 1."
            prm.write("{} ! {}\n".format(self.output_sampling, string_10))
            string_11 = "Image frame offset in pixels of the input wave. The parameter is used " \
                        "only if the Flag in line 09 is set to 1."
            prm.write(str(self.img_frame_offset[0]) + ', ' + str(self.img_frame_offset[1]) +
                      "{}, {} ! {}\n".format(self.img_frame_offset[0], self.img_frame_offset[1],
                                             string_11))
            string_12 = "Image frame rotation in [deg] with respect to the input wave " \
                        "horizontal axis. The parameter is used only if the Flag in line 09 is " \
                        "set to 1."
            prm.write("{} ! {}\n".format(self.img_rot, string_12))
            string_13 = "Coherence calculation model switch: 1 = averaging of coherent sub " \
                        "images explicit focal variation but quasi-coherent spatial envelope, " \
                        "2 = averaging of coherent sub images with explicit focal and angular " \
                        "variation, 3 = quasi-coherent linear envelopes, 4 = Fourier-space " \
                        "synthesis with  partially coherent TCC, 5: averaging of coherent sub " \
                        "images with explicit  focal, angular, and frozen lattice variation)."
            prm.write("{} ! {}\n".format(self.coherence_model, string_13))
            string_14 = "Flag and parameters for partial temporal coherence: <ptcflg> = flag (" \
                        "0=OFF, 1=ON), <f-spread> = focus spread (1/e) half width [nm]."
            prm.write("{}, {} ! {}\n".format(self.temp_coherence[0], self.temp_coherence[1],
                                             string_14))
            string_15 = "Flag and parameters for partial spatial coherence: <pscflg> = flag (" \
                        "0=OFF, 1=ON), <s-conv> = beam convergence (1/e) half width [mrad]."
            prm.write("{}, {} ! {}\n".format(self.spat_coherence[0], self.spat_coherence[1],
                                             string_15))
            string_16 = "Flag and parameters for applying the detector MTF: <mtfflag> = flag (" \
                        "0=OFF, 1=ON), <mtf-scale> = calculation scale of the mtf = (sampling " \
                        "rate experiment)/(sampling rate simulation), <mtf-file> = File name " \
                        "string to locate the MTF data. Use quotation marks when the string " \
                        "includes space characters."
            prm.write("{}, {}, '{}' ! {}\n".format(self.mtf[0], self.mtf[1], self.mtf[2],
                                                   string_16))
            string_17 = "Flag and parameters for a vibration envelope: <vibflg> = flag (0=OFF, " \
                        "1=ON-ISO, 2=ON-ANISO), <vibprm1>, <vibprm1> = vibration RMS amplitudes " \
                        "[nm], <vibprm3> = orientation [deg] of the primary vibration amplitude " \
                        "w.r.t. the horizontal image axis."
            prm.write("{}, {}, {}, {} ! {}\n".format(self.vibration[0], self.vibration[1],
                                                     self.vibration[2], self.vibration[3],
                                                     string_17))
            string_18 = "Number of aberration definitions following this line."
            prm.write("{} ! {}\n".format(self.number_of_aberrations, string_18))
            for key in self.aberrations_dict:
                prm.write('{} {:.4f} {:.4f} ! {}\n'.format(key, self.aberrations_dict[key][0],
                                                           self.aberrations_dict[key][1],
                                                           aberrations[key]))
            #prm.write("{} ! {}\n".format(self.number_of_aberrations, string_18))
            #for key in self.aberrations_dict:
            #    prm.write('{}, {}, {} ! {}\n'.format(key, self.aberrations_dict[key][0],
            #                          self.aberrations_dict[key][1], aberrations[key]))
            string_19 = "Objective aperture radius [mrad]. Set to very large values to deactivate."
            prm.write("{} ! {}\n".format(self.oa_radius, string_19))
            string_20 = "Center of the objective aperture with respect to the zero beam [mrad]."
            prm.write("{}, {} ! {}\n".format(self.oa_position[0], self.oa_position[1], string_20))
            string_21 = "Number variable of loop definitions following below."
            prm.write("{} ! {}".format(self.number_of_loops, string_21))

        # Sort prm file
        with open(prm_filename, 'r+') as prm:
            content = prm.readlines()

            length = []
            for line in content:
                length.append(len(line.rsplit('!', 1)[0]))

            spacer = np.max(length)
            prm.seek(0)
            prm.truncate()
            for line in content:
                line = line.split('!')[0].ljust(spacer) + ' !' + line.split('!')[1]
                prm.write(line)
        if output:
            print("Parameters successfully saved to file '{}'!".format(prm_filename))
