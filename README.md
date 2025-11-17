# pyCallisto

pyCallisto is an image processing and data analysis Python library developed for e-CALLISTO data (http://www.e-callisto.org/). This library can be used to process data obtained using various solar radio spectrometers available around the globe.

# Installation

You can install pyCallisto using pip:

```bash
pip install pyCallisto
```

For more details, visit: https://pypi.org/project/pyCallisto/

# Prerequisites and Standard Python Libraries 

**Required:**
- Python 3.6 or higher
- numpy
- matplotlib
- astropy

**Optional (for GOES integration):**
- sunpy
- scipy

Install optional dependencies with:
```bash
pip install sunpy scipy
```

# Features

- **Spectrogram Visualization**: Generate spectrograms from FITS files with customizable colormaps and parameters
- **Time and Frequency Slicing**: Extract specific time ranges and frequency bands from observations
- **Light Curve Analysis**: Create mean light curves by collapsing data along time axis
- **Spectrum Analysis**: Generate frequency spectra for specific time instances
- **GOES Integration**: Overlay GOES X-ray flux data on spectrograms for comprehensive solar event analysis
- **Background Subtraction**: Automated background estimation and removal
- **Universal Plot**: Combined visualization of spectrogram, light curve, and spectrum
- **File Concatenation**: Merge multiple FITS files along time axis for extended observations

# Contributors

**Mr. Ravindra Pawase**  
Data Scientist  
Cummins Inc.  
Pune-411 004, India

**Dr. K. Sasikumar Raja**  
Assistant Professor  
Indian Institute of Astrophysics  
II Block, Koramangala, Bengaluru-560 034, India

**Mr. Mrutyunjaya Muduli**  
B.E. Computer Science and Engineering  
Department of Computer Science and Engineering  
HKBK College of Engineering  
22/1, Opposite Manyata Tech Park, Nagavara, Bengaluru-560 045, India

# Feedback

For feedback, queries, or feature requests, please contact:

- ravi.pawase@gmail.com
- sasikumarraja@gmail.com
- mudulimrutyunjaya42@gmail.com

# Citation

If you find pyCallisto useful in your work, we appreciate acknowledgment. We recommend using the following citation:

"This work makes use of the pyCallisto library, which is available at https://github.com/ravipawase/pyCallisto"

