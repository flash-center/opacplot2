# opacplot2

[![Build Status](https://travis-ci.org/flash-center/opacplot2.svg?branch=master)](https://travis-ci.org/flash-center/opacplot2)

Python package for manipulating Equation of State (EoS) and Opacity data.

`Opacplot2` comes with an EoS Table conversion tool named [`opac-convert`](#opac-convert).
It also comes with an EoS Table comparison tool named [`opac-error`](#opac-error).

### Dependencies

`opacplot2`'s dependencies include:

* numpy 
* nose 
* six 
* pytables 
* matplotlib 
* scipy
* hedp (https://github.com/luli/hedp)

They can be installed as follows:

```shell
pip install numpy nose six pytables matplotlib scipy
pip install git+https://github.com/luli/hedp
```

### Installation 

This module requires Python 2.7 or 3.5. The latest version can be installed with

```shell
pip install git+https://github.com/flash-center/opacplot2
```

If you have the Propaceos Python reader, in order to include it in the 
installation, you must install `opacplot2` as follows:

```shell
git clone https://github.com/flash-center/opacplot2
cp /path/to/opg_propaceos.py opacplot2/opacplot2/
cd opacplot2
python setup.py install
```

### Documentation

Full documentation can be found in the `doc/` directory (html or text).

---

<a name="opac-convert"></a>
# opac-convert

Command line tool for converting EoS Table formats into the IONMIX format
that comes with `opacplot2`.

Supported input file formats:

* Propaceos (not distributed, contact jtlaune at uchicago dot edu.)
* SESAME (.ses)
* QEOS SESAME (.mexport)
* MULTI (.opp, .opr, .opz, .eps)

Currently the only supported output format is IONMIX.

### Usage

```bash
opac-convert [options] myfile.ext
```

`opac-convert` will attempt to read your file extension and convert it to
IONMIX accordingly.
If it is unable to read the extension, you can use the input flag to specify
your filetype.
Some files need additional information to write to IONMIX,
such as atomic numbers. These you must specify with the command line options
shown below.

### Options

| Option | Action |
|:-------|--------|
|-i, --input| Specify the input filetype (`propaceos`, `sesame`, `multi`)|
|--Znum| Comma separated list of atomic numbers.|
|--Xfracs| Comma separated list of element fractions.|
|--outname| Specify the output filename.|
|--log| Comma separated list of logarithmic data.|
|--tabnum| SESAME table number (defaults to last).|

### Example

To specify the files atomic numbers, one may use `--Znum` with a comma separated
list of integers. If more than one atomic number is given, 
one must also specify the element fractions with `--Xfracs`.
For example, take a SESAME table for CH named `myfile.ses`:

```bash
opac-convert --Znum 1, 6 --Xfracs .5, .5 myfile.ses
```

This will convert `myfile.ses` to an IONMIX file named `myfile.cn4`.

### Logarithmic Data 

If you would like to take the log of the data before you write it to the IONMIX
file, use `--log` with a comma separated list of the data keys as shown below.
Each key specified will be written to IONMIX after the base 10 logarithm has
been applied.

| Data | key |
|------|-----|
|Ion number density|`idens`|
|Temperatures|`temps`|
|Average ionization|`Zf_DT`|
|Ion pressure|`Pi_DT`|
|Electron pressure|`Pec_DT`|
|Ion internal energy|`Ui_DT`|
|Electron internal energy|`Uec_DT`|
|Opacity bounds|`groups`|
|Rosseland mean opacity|`opr_mg`|
|absorption Planck mean opacity|`opp_mg`|
|emission Planck mean opacity|`emp_emg`|


For example, in order to specify that the emission Planck Mean Opacity be written
logarithmically:

```bash
opac-convert --log emp_mg my-file.ext
```

### Troubleshooting

#### Invalid Literal for `int()`

The `--log` flag may be used to fix the following error: 

```bash
ValueError: invalid literal for int() with base 10
```

This error arises when the exponent for the data is more than 2 digits long, 
which IONMIX does not support. 
What that usually means is that the data was originally
stored logarithmically and must be written back to IONMIX as logarithmic data.

---

<a name="opac-error"></a>
# opac-error

Command line tool for comparing two EoS/Opacity tables.

Supported input file formats:

* Propaceos (not distributed, contact jtlaune at uchicago dot edu.)
* SESAME (.ses)
* QEOS SESAME (.mexport)
* IONMIX

### Usage

```bash
opac-convert [options] myfile_1.ext myfile_2.ext
```

Much like `opac-convert`, this tool will first attempt to read the extensions from your EoS tables
in order to open them up. However, some file types require additional information, as in 
`opac-convert`, which can be specified in the `[options]`. The options have many of the same
names as `opac-convert` but suffixed by `_1` for file 1 or `_2` for file 2.

Then, `opac-error` will create error reports for the following data:

* Average ionization
* Electron Ppessure
* Ion pressure
* Electron energy
* Ion energy

Included in the error report, we have

* Root mean squared % error
* Maximum absolute % error

If the `--plot` flag is called, `opac-convert` will also make error plots and save them as images
to the current directory.

### Options

| Option | Action |
|:-------|--------|
|--filetypes| Comma separated list of file types.|
|--mpi_#| Mass per ion (g) for file 1 or 2.|
|--Znum_#| Comma separated list of atomic numbers for file 1 or 2.|
|--Xfrac_#| Comma separated list of number fractions for file 1 or 2.|
|--filters_#| `dens_filter, temp_filter` for file 1 or 2.|
|--tabnum_#| SESAME table number for file 1 or 2.|
|--plot| Create % error plots for data.|
|--writelog| Write log file with % errors for data.|
|--lin_grid| Plot using linear axes.|

### How It Works

First, `opac-convert` takes an intersection of the two dens/temp grids from each file.
Then it interpolates the data from both files onto the intersection dens/temp grids and is then
able to create an error report.