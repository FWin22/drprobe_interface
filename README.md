# drprobe_interface

drprobe_interface provides a Python interface for using the command line tools of the [Dr Probe
software
package.](http://www.er-c.org/barthel/drprobe/)

## Requirements
[Dr Probe Command Line Tools](http://www.er-c.org/barthel/drprobe/drprobe-download.html) must be installed.

## Example
Import the DrProbe package:
    
    import drprobe as drp

Run the celslc tool:
    
    drp.commands.celslc(cel_file, slice_name, ht, *kwargs)
    
Create a msa-parameter file for the multislice algorithm:

    msa_prm = drp.msaprm.MsaPrm()
    
You can edit the parameters of the file, for example the object tilt:

    msa_prm.tilt_x = 1    # Applies 1° tilt in x-direction
    msa_prm.tilt_y = 0.2  # Applies 0.2° tilt in y-direction
    
Then save the parameter file:

    msa_prm.save_msa_prm('msa.prm')
    
Alternatively, you can also load an existing parameter file and edit it:

    msa_prm = dr.msaprm.load_msa_prm('msa.prm')
    
Execute the multislice algorithm:

    drp.commands.msa('msa.prm', 'output_name', *kwargs)

Create (or load) a wavimg-parameter file and edit the parameters accordingly, for example setting the output dimensions to 256 x 256 pixels. Afterwards, save the parameter file:

    wavimg_prm = drp.wavimgprm.WavimgPrm()
    wavimg_prm.output_dim = (256, 256)
    wavimg_prm.save_wavimg_prm('wavimg.prm')
    
Simulate the image:

    drp.commands.wavimg('wavimg.prm', *kwargs)
    
  
