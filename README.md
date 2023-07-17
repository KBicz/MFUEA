# MFUEA
Program to model the modulated stellar flare by the stellar rotation using the Differential Evolution method (DE). To properly use the program you have to install gnuplot on your computer. Compile it using the command "cargo build --release". The input file to the program needs the four-column file where:
- the first column is the time,
- the second column is the flux of the star with rotational modulation,
- the third column is the error of the flux,
- the fourth column is the model of rotational modulation,

Additionally, the first line of the input file needs to start with the comment: # index_of_the_flare_start_in_the_data index_of_the_flare_end_in_the_data (example: # 81 141). An example of how to use the software is in the "example" directory.

    Program mfuea written by K. Bicz for Linux, Mac, and Windows, version of 13 Mar 2023.
    Usage: mfuea <-lc=file> [-cross=f64] [-mut=f64] [-npop=int] [-niter=f64] [-prot=f64]
                 [-imin=f64] [-imax=f64] [-lamin=f64] [-lamax=f64] [-lomin=f64] [-lomax=f64]
                 [-nparams=4/7/8] [-npts=i64] [-teff=f64] [-rad=f64] [-terr=f64] [-raderr=f64]
                 [-logg=f64] [-respf=str] [-frad=f64] [-famp] [--limbd] [--noprof] [--energy]
                 [--plot] [--save]
 
         option  -lc      : file with the light curve of the event.
                 -cross   : crossing value (default cross = 0.1).
                 -mut     : mutation value (default mut = 0.0).
                 -npop    : number of population (default npop = (nparams+3)*10).
                 -niter   : number of iterations with mutation and crossing (default niter = 10000).
                 -imin    : minimal inclination angle of the star (default imin = 10).
                 -imax    : maximal inclination angle of the star (default imax = 90).
                 -lamin   : minimal latitude of the flaring region (default lamin = -90).
                 -lamax   : maximal latitude of the flaring region (default lamax = 90).
                 -lamin   : minimal longitude of the flaring region (default lamin = 0, east limb).
                 -lamax   : maximal longitude of the flaring region (default lamax = 360).
                 -nparams : number of parameters to fit profile (default nparams = 4).
                 -npts    : make region not the point but area made of npts points (default npts = 100)
                 -prot    : rotational period of the star (default prot = 0.2353 days).
                 -teff    : effective temperature of the star (default teff = 4885.0 K).
                 -terr    : error of the effective temperature of the star (default terr = 128.079 K
                 -rad     : radius of the star (default rad = 0.799503 Rsun
                 -raderr  : error of the radius of the star (default rad = 0.0540905 Rsun
                 -logg    : log(g) of the star (def logg = 4.53552).
                 -frad    : radius of the flaring region in degrees (by default it is estimated).
                 -famp    : force the new amplitude of the flare (by default it is taken from data).
                 -respf   : file with the response function of the observational instrument
                 --limbd  : use limb darkening in the fit.
                 --noprof : assume that fitted profile from lc file is the final profile of the event.
                 --energy : estimate the energy of the flare before and after the correction.
                 --plot   : display the results.
                 --save   : save the results to file.
