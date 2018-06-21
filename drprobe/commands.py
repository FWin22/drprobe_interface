# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 09:56:58 2016

@author: fwinkler
"""

import os
import re
import subprocess


def cellmuncher(cel_file, output_file, attach_cel=None, attach_direction=None,
                repeat=None, sort=None, cif=False, override=False):
    """
    Runs cellmuncher. Supports only a few basic options at the moment.

    Parameters
    ----------
    cel_file : str
        Super-cell file
    output_file : str
        New cel file
    attach_cel : str, optional
        cel file of the cel that will be attached to cel_file
    attach_direction : str, optional
        direction in which the cel shall be attached (x, y, z)
    repeat : (str, int), optional
        repeat the cel file in a certain direction, e.g., ('x', 5) repeats the cel file 5 times in
        x-direction)
    sort : list of str, optional
        sort cel file according to strings (eg, 'x' for x-coordinate or 'e' for element)
    cif : bool, optional
        Export cif file
    override : bool, optional
        Overrides output file without asking
    """

    _cellmuncher_options = {}

    if cif:
        _cellmuncher_options['cif'] = ' --cif'
    if attach_cel:
        _cellmuncher_options['attach_cel'] = ' --attach-cell={},XMS,{}'.format(attach_cel,
                                                                               attach_direction)
    if repeat:
        _cellmuncher_options['repeat'] = ' --repeat={},{}'.format(repeat[0], repeat[1])

    if sort:
        for item in sort:
            _cellmuncher_options['sort'] += ' -s={}'.format(item)
    if override:
        _cellmuncher_options['override'] = ' --override'

    # Define command
    command = "cellmuncher -f {} -o {}".format(cel_file, output_file)

    for key in _cellmuncher_options:
        command += _cellmuncher_options[key]

    # Run the celslc command
    subprocess.call(command, shell=True)

    print('Performed cellmuncher with the following command:\n', command)


def celslc(cel_file, slice_name, ht, nx=None, ny=None, nz=None, abf=None, absorb=False,
           dwf=False, buni=None, fl=False, nv=None, pot=False, pps=False, prj=None, rev=False,
           ssc=None, tla=None, _3dp=False, inf=None, rti=False, silent=False, output=False):
    """
    Runs celslc from Dr. Probe.
    Requires installation of Dr Probe command line tools.
    For further information please read the 'celslc howto.txt' file from Dr. probe.

    Parameters
    ----------
    cel_file : str
        CEL structure file file (e.g. "atoms.cel")
        Can also be a cif structure file as defined in [http://www.iucr.org/resources/cif].
    slice_name : str
        Slice file name (e.g. "slice")
        The format of the output slice files is uncompressed binary raw data.
    ht : int
        Electron energy in keV (10 ... 1300)
    nx : int, optional
        Vertical cell sampling (32 ... 2048)
        Not required when loading external 3D potentials. nx must be defined in the potential
        header.
    ny : int, optional
        Horzontal cell sampling (32 ... 2048)
        Not required when loading external 3D potentials. ny must be defined in the potential
        header.
    nz : int, optional
        Number of slices (1 ... 2048)
        Not required when loading external 3D potentials. nz must be defined in the potential
        header.
    abf : float, optional
        Apply absorptive form factors as a fraction of the elastic form factors. The fraction
        is specified by the input number.
        The option abf can only be used in combination with the option dwf.
    absorb : bool, optional
        Apply built-in absorptive form factors calculated according to A. Weickenmeier and H.
        Kohl [Acta Cryst. A47 (1991) 590-597].
        The option absorb can only be used in combination with the option dwf.
        Using the option absorb will override and deactivate the option abf.
        absorb is identical to the -abs command in Dr. Probe.
    dwf : bool, optional
        Apply Debye-Waller factors as specified in the cel_file.
        The electron scattering potentials are convoluted by Debye-Waller factors.
    buni : float, optional
        Specify uniform Debye-Waller factor.
        Resets all B_ISO Parameters to the given value.
        The option buni can only be used in combination with dwf.
    fl : bool, optional
        Frozen lattice simulation.
        For frozen-lattice simulations you should combine this option with the option nv.
        The option fl will override and deactivate the option dwf.
    nv : int, optional
        Number of frozen lattice variants per slice (1 ... 2048)
    pot : bool, optional
        Export projected electron scattering potentials for each slice to files *.pot.
        Potentials are stored without the relativistic correction.
    pps : bool, optional
        Output projected potentials instead of phase gratings to the slice files.
        Potentials are stored without the relativistic correction.
    prj : tuple, optional
        Super-cell re-orientation (u, v , w, u, v, w, a, b, c).
        Further information in the 'celcls howto.txt'.
    rev : bool,  optional
        Slicing in reversed sequence
    ssc : int, optional
        This option allows to calculate the phase grating for one specific slice of the object.
        The slice is identified by an index number. The indexing starts form number 1 and goes
        up to the number of slices defined with the nz option.
    tla : tuple, optional
        Additional shift of atoms (x, y, z).
        This option can be used to shift the whole structure in the chosen final orthorhombic
        super-cell by a given amount in fractional coordinates. All atoms are shifted by the
        same amount and their final fractional coordinates are wrapped back periodically into
        the range [0, 1).
    _3dp : bool, optional
        Activates 3D potential simulation.
    inf : int, optional
        Read external 3D potentials.
        By setting inf to 10 you select the 4th column as the data columns. Setting inf to 11
        selects the 5th column, and so on.
    rti : bool, optional
        Run time information.
    silent : bool, optional
        If True, fully deactivate terminal output.
    output : bool, optional
        Activates terminal output of celslc command.
    """

    _celslc_options = {}

    # If necessary, create folder for output slices.
    directory = os.path.split(slice_name)[0]
    if directory:
        if not os.path.isdir(directory):
            os.makedirs(directory)

    if rev:
        _celslc_options['rev'] = ' -rev'
    if fl:
        _celslc_options['fl'] = ' -fl'
    if nv:
        _celslc_options['nv'] = ' -nv {}'.format(nv)
    if dwf:
        _celslc_options['dwf'] = ' -dwf'
    if buni:
        _celslc_options['buni'] = ' -buni {}'.format(buni)
    if absorb:
        _celslc_options['absorb'] = ' -abs'
    if abf:
        _celslc_options['abf'] = ' -abf {}'.format(abf)
    if pot:
        _celslc_options['pot'] = ' -pot'
    if _3dp:
        _celslc_options['_3dp'] = ' -3dp'
    if inf:
        _celslc_options['inf'] = ' -inf {}'.format(inf)
    if pps:
        _celslc_options['pps'] = ' -pps'
    if ssc:
        _celslc_options['ssc'] = ' -ssc {}'.format(ssc)
    if rti:
        _celslc_options['rti'] = ' -rti'
    if silent:
        _celslc_options['silent'] = ' -silent'
    if prj:
        _celslc_options['prj'] = ' -prj'
        for i in prj:
            _celslc_options['prj'] += ' {},'.format(i)
        _celslc_options['prj'] = _celslc_options['prj'][:-1]
    if tla:
        _celslc_options['tla'] = ' -tla'
        for i in tla:
            _celslc_options['tla'] += ' {},'.format(i)
        _celslc_options['tla'] = _celslc_options['tla'][:-1]

    # Define command
    if cel_file.endswith('.cel') or cel_file.endswith('.txt'):
        command = "celslc -cel {} ".format(cel_file)
    elif cel_file.endswith('.cif'):
        command = "celslc -cif {} ".format(cel_file)
    else:
        command = "celslc -cel {} ".format(cel_file)

    if inf is None:
        command += "-slc {} -nx {} -ny {} -nz {} -ht {}".format(slice_name, nx, ny, nz, ht)
    else:
        command += "-slc {} -ht {}".format(slice_name, ht)

    for key in _celslc_options:
        command += _celslc_options[key]

    # Run the celslc command
    subprocess.call(command, shell=True)

    if output:
        print('Performed celslc with the following command:\n', command)


def msa(prm_file, output_file, input_image=None, inw=None, px=None, py=None, lx=None, ly=None,
        foc=None, tx=None, ty=None, otx=None, oty=None, sr=None, abf=None, buni=None, uuni=None,
        ctem=False, txtout=False, _3dout=False, gaussap=False, wave=False, avwave=False,
        detimg=False, verbose=False, debug=False, lapro=False, waveft=False, avwaveft=False,
        vtx=None, silent=False, rti=False, output=False):
    """
    Runs msa from Dr. Probe

    Parameters
    ----------
    prm_file : str
        The filename of the msa parameter file.
    output_file : str
        The filename of the output data. Example: 'Test.wav'.
        The result is a 4-byte floating point array for a STEM image calculation and an 8-byte
        complex array for a CTEM wave calculation.
    input_image : str, optional
        Input image filename.
        This option is used for the application of partial spatial coherence in STEM mode
        calculations only.
        Equivalent to the option '-in' in msa in Dr. Probe.
    inw : (str, int), optional
        Tuple. (Input wave filename, Input slice number).
        Equivalent to the option '-inw' in msa in Dr. Probe.
    px : int, optional
        Horizontal scan pixel number.
        Defines the x-scan position or scan column number starting with 0 up to number of scan
        columns -1.
    py : int, optional
        Vertical scan pixel number.
        Defines the y-scan position or scan row number starting with 0 up to number of scan rows
        -1.
    lx : int, optional
        Last horizontal scan pixel.
        Defines the last x-scan position (scan column) number starting from 0 up to # scan
        columns -1 to be calculated, in this case -px defines the first scan column to be
        calculated.
    ly : int, optional
        Last vertical scan pixel.
        Defines the last y-scan position (scan row) number starting from 0 up to # scan
        rows -1 to be calculated, in this case -py defines the first scan row to be calculated.
    foc : float, optional
        Defines the probe or image wave defocus in nm. Overrides the value set in the parameter
        file.
    tx : float, optional
        Defines the x beam tilt in mrad. Overrides the value set in the parameter file.
    ty : float, optional
        Defines the y beam tilt in mrad. Overrides the value set in the parameter file.
    otx : float, optional
        Defines the x object tilt in mrad. Overrides the value set in the parameter file.
    oty : float, optional
        Defines the y object tilt in mrad. Overrides the value set in the parameter file.
    sr : float, optional
        Effective source radius in nm.
        Defines the half-width of the effective geometrical source profile in nm. Overrides the
        value set in the parameter file. Half-width definitions depend on the source profile
        chosen in the parameter file.
    abf : float, optional
        Applies a fractional absorption potential. Only when the slice files are given with
        potential data, not possible in case of phase grating data.
    buni : float, optional
        Applies a universal Debye-Waller factor. The float value is the isotropic Debye-Waller
        parameter (Biso) in nm^2. Only when the slice files are given with potential data,
        not possible in case of phase grating data.
    uuni : float, optional
        Applies a universal Debye-Waller factor. The float value is the isotropic Debye-Waller
        parameter (Uiso) in A^2. Only when the slice files are given with potential data,
        not possible in case of phase grating data.
    ctem : bool, optional
        Activates CTEM mode calculating a coherent(!) wave function. Image aberrations can be
        applied to the exit-plane wave function.
    txtout : bool, optional
        STEM mode only. Outputs the simulation result in form of a text file, which includes a
        header and lists of data. See 'msa howto.txt' (section 4) for an explanation of the text
        file format.
    _3dout : bool, optional:
        STEM mode only. The output of results result is in 4-byte floating point arrays. The
        third dimension represents the thickness series. One file is generated for each detector
        and a respective suffix with the detector name is appended. This option applies also
        when using input image files with the option 'input_image' ('-in' in msa). Files
        will be handled consistently only, i.e. 3d -> 3d or 2d -> 2d.
    gaussap : bool, optional
        Applies a gaussian illumination aperture in STEM mode leading to a larger delocalization
        of the incident probe. The 1/e width of the applied Gaussian is defined by the aperture
        size parameter in mrad.
    wave : bool, optional
        Activates wave export at selected slices. Wavefunctions are saved to file name = name of
        output file. The file name is extended with ".wav". Index suffixes are added to denote
        the scan pixel, the object slice and the frozen phonon variation index.
    avwave : bool, optional
        Activates the output of an average wavefunctions over multiple frozen lattice calculations.
    detimg : bool, optional
        Images representing the applied detector functions in diffraction space are saved to files.
    verbose : bool, optional
        Activates additional program output to the console.
    debug : bool, optional
        Activates additional debug output to the console.
    lapro : bool, optional
        Use large-angle propagators exp(-2*Pi/Lambda*(1/cos(theta)-1)) instead of the default
        Fresnel propagators exp(-I*Pi*Lambda*q^2).
    waveft : bool, optional
        Skip the extra inverse Fourier transform done when using option "wave".
    avwaveft : bool, optional
        Skip the extra inverse Fourier transform done when using option "avwave"
    vtx : int, optional
        Enable the use of vortex probes in STEM mode. The integer specifies the orbital angular
        momentum of the probe.
    silent : bool, optional
        Flag for deactivating all console output.
    rti : bool, optional
        Run time information.
    output : bool, optional
        Flag for terminal output
    """

    _msa_options = {}

    # Make folder for the output files if it doesn't exist already
    directory = os.path.split(output_file)[0]
    if directory:
        if not os.path.isdir(directory):
            os.makedirs(directory)

    if input_image:
        _msa_options['input_image'] = ' -in {}'.format(input_image)
    if inw:
        _msa_options['inw'] = ' -inw {} {}'.format(inw[0], inw[1])
    if px:
        _msa_options['px'] = ' -px {}'.format(px)
    if py:
        _msa_options['py'] = ' -py {}'.format(py)
    if lx:
        _msa_options['lx'] = ' -lx {}'.format(px)
    if ly:
        _msa_options['ly'] = ' -ly {}'.format(py)
    if foc:
        _msa_options['foc'] = ' -foc {}'.format(foc)
    if tx:
        _msa_options['tx'] = ' -tx {}'.format(tx)
    if ty:
        _msa_options['ty'] = ' -ty {}'.format(ty)
    if otx:
        _msa_options['otx'] = ' -otx {}'.format(otx)
    if oty:
        _msa_options['oty'] = ' -oty {}'.format(oty)
    if sr:
        _msa_options['sr'] = ' -sr {}'.format(sr)
    if abf:
        _msa_options['abf'] = ' -abf {}'.format(abf)
    if buni:
        _msa_options['buni'] = ' -buni {}'.format(buni)
    if uuni:
        _msa_options['uuni'] = ' -uuni {}'.format(uuni)
    if ctem:
        _msa_options['ctem'] = ' /ctem'
    if txtout:
        _msa_options['txtout'] = ' /txtout'
    if _3dout:
        _msa_options['_3dout'] = ' /_3dout'
    if gaussap:
        _msa_options['gaussap'] = ' /gaussap'
    if wave:
        _msa_options['wave'] = ' /wave'
    if avwave:
        _msa_options['avwave'] = ' /avwave'
    if detimg:
        _msa_options['detimg'] = ' /detimg'
    if verbose:
        _msa_options['verbose'] = ' /verbose'
    if debug:
        _msa_options['debug'] = ' /debug'
    if lapro:
        _msa_options['lapro'] = ' /lapro'
    if waveft:
        _msa_options['waveft'] = ' /waveft'
    if avwaveft:
        _msa_options['avwaveft'] = ' /avwaveft'
    if vtx:
        _msa_options['vtx'] = ' /vtx {}'.format(vtx)
    if silent:
        _msa_options['silent'] = ' /silent'
    if rti:
        _msa_options['rti'] = ' /rti'

    # Define command
    command = "msa -prm {} -out {}".format(prm_file, output_file)
    for key in _msa_options:
        command += _msa_options[key]

    # Run msa command
    subprocess.call(command, shell=True)

    if output:
        print('Performed msa with the following command:\n', command)


def wavimg(prm_file, output_file=None, foc=None, btx=None, bty=None, oar=None,
           sbshx=None, sbshy=None, sil=False, dbg=False, nli=False, rnsb=False,
           rti=False, output=False):
    """
    Runs wavimg from Dr. Probe

    Parameters
    ----------
    prm_file : str
        The filename of the wavimg parameter file.
    output_file : str, optional
        The filename of the output data.
    foc : float, optional
        Image defocus in nm.
    btx : float, optional
        Beam tilt x in mrad.
    bty : float, optional
        Beam tilt y in mrad.
    oar : float, optional
        Objective aperture radius in mrad.
    sbshx : float, optional
        Apply MTF for a shifted sideband in x, sbshx is given in pixels
    sbshy : float, optional
        Apply MTF for a shifted sideband in y, sbshy is given in pixels
    sil : bool, optional
        Deactivates console output.
    dbg : bool, optional
        Activates extra debug consolue output
    nli : bool, optional
        no loop image output in map output mode
    rnsb : bool, optional
        Using to renorm the MTF after sideband shift.
    rti : bool, optional
        Run time information.
    output : bool, optional
        Flag for terminal output
    """

    _wavimg_options = {}

    # Check if output_file is given as parameter
    if output_file:
        directory = os.path.split(output_file)[0]
        _wavimg_options['output_file'] = ' -out {}'.format(output_file)
    else:
        with open(prm_file, 'r') as prm:
            _content = prm.readlines()
            _content = [re.split(r'[,\s]\s*', line) for line in _content]
            directory = os.path.split(_content[5][0])[0].replace("'", "")

    if btx:
        _wavimg_options['btx'] = ' -btx {}'.format(btx)
    if bty:
        _wavimg_options['bty'] = ' -bty {}'.format(bty)
    if foc:
        _wavimg_options['foc'] = ' -foc {}'.format(foc)
    if oar:
        _wavimg_options['oar'] = ' -oar {}'.format(oar)
    if sbshx:
        _wavimg_options['sbshx'] = ' -sbshx {}'.format(sbshx)
    if sbshy:
        _wavimg_options['sbshy'] = ' -sbshy {}'.format(sbshy)
    if sil:
        _wavimg_options['sil'] = ' /sil'
    if dbg:
        _wavimg_options['dbg'] = ' /dbg'
    if nli:
        _wavimg_options['nli'] = ' /nli'
    if rnsb:
        _wavimg_options['rnsb'] = ' /rnsb'
    if rti:
        _wavimg_options['rti'] = ' /rti'

    # Make folder for outputfiles if it doesn't exist already
    if directory:
        if not os.path.isdir(directory):
            os.makedirs(directory)

    # Define command
    command = "wavimg -prm {}".format(prm_file)

    for key in _wavimg_options:
        command += _wavimg_options[key]

    # Run wavimg command
    subprocess.call(command, shell=True)

    if output:
        print('Performed wavimg with the following command:\n', command)
