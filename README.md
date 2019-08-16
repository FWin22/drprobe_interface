# drprobe_interface

drprobe_interface provides a Python interface for using the command line tools of the [Dr Probe
software
package.](http://www.er-c.org/barthel/drprobe/)

## Requirements
[Dr Probe Command Line Tools](http://www.er-c.org/barthel/drprobe/drprobe-download.html) must be installed.

## Installation

    pip install drprobe-interface

## Example
The user is referred to the documentation of the Dr Probe command line tools for further
information on S(TEM) image simulation. The following example provides a basic overview on the
workflow of using the drprobe_interface package.
A more detailled example for an HRTEM image simulation of SrTiO3 is given in the [examples folder](https://github.com/FWin22/drprobe_interface/tree/master/examples).

Import the drprobe_interface package:
    
    import drprobe as drp

Run the celslc tool:
    
    drp.commands.celslc(cel_file, slice_name, ht, *kwargs)
    
Create a msa-parameter file for the multislice algorithm:

    msa_prm = drp.msaprm.MsaPrm()

Alternatively, you can also load an existing parameter file and edit it:

    msa_prm = dr.msaprm.load_msa_prm('msa.prm')
    
You can edit the parameters of the file, for example the object tilt:

    msa_prm.tilt_x = 1    # Applies 1° tilt in x-direction
    msa_prm.tilt_y = 0.2  # Applies 0.2° tilt in y-direction
    
Then save the parameter file:

    msa_prm.save_msa_prm('msa.prm')
    
Execute the multislice algorithm:

    drp.commands.msa('msa.prm', 'output_name', *kwargs)

Create (or load) a wavimg-parameter file and edit the parameters accordingly, for example setting
the output dimensions to 256 x 256 pixels. Afterwards, save the parameter file:

    wavimg_prm = drp.wavimgprm.WavimgPrm()
    wavimg_prm.output_dim = (256, 256)
    wavimg_prm.output_files = 'img.dat'
    wavimg_prm.save_wavimg_prm('wavimg.prm')
    
Simulate the image:

    drp.commands.wavimg('wavimg.prm', *kwargs)

Depending on the datatype and size, the simulated image can be loaded in python using numpy:

    import numpy as np
    img = np.fromfile('img.dat', dtype='float32')
    # Change dtype to 'complex64', when simulating complex-valued wavefunctions
    img = img.reshape(wavimg_prm.output_dim)
  
