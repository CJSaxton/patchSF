# patchSF
Given piecewise patchy time series data, calculate the Kolmogorov structure functions at given order, and respective uncertainties.  The code "patchsf" calculates the Kolmogorov structure function (at given order) from discrete time series data fron a 4-column input file.  Write a 4-column output file containing timescales, SF denominator, SF numerator, uncertainty.

The first command-line argument is the input data file.

INPUT FORMAT:
The input file has four columns:

    t       = time of observation (temporal bin centroid)
    dt      = width of temporal bin (or time uncertainty)
    z       = signal or count-rate during this time bin
    dz      = uncertainty in the signal or count-rate
    
OUTPUT FORMAT:

    tau     = timescale of SF evaluation
    W       = a weight or normalisation denominator in the SF
    W S     = numerator in the structure function at timescale tau
    W dS    = uncertainty propagated into the product (W S)
    
Thus the structure function at timescale tau, and its uncertainty are:

    S_n(tau)  = (W S) / W   =  (column 3) / (column 2)
    dS_n(tau) = (W dS) / W  = (column 4) / (column 2)

The code writes the denominator "W" as a separate column as a precaution.  With temporally patchy data, in some numerical methods, there are sometimes circumstances when W=0 at some timescale tau.  Keeping "W" distinct enables an informed user to discard those rows from the output file.

The sequence of tau values is the union of a linearly spaced sequence (tmin<=tau<=tmax) and a logarithmically spaced sequence in the same domain.  To increase or decrease the resolution of this output grid, vary the "try" command-line parameter.


EXAMPLE:

    ./patchSF signal.dat start=1. end=20. n=2 tmin=10. tmax=1.e3 ^results.dat
 
 COMMAND-LINE PARAMETERS & OPTIONS:
 
    start   = use input data after this time
    end     = use input data before this time
    -logx   flag: take logarithm of independent variable (time)
    -logy   flag: take logarithm of dependent variable (signal)
    
    n       = order of the structure function
    try     = number of tau grid points -1
    tmin    = shortest time-scale (tmin <= tau)
    tmax    = longest time-scale (tau <= tmax)
    
Alternatively, calculate a Leahy normalised power spectrum.
    -pds    flag: compute power density spectrum, not structure function
 
 Results are written into a file named after the host machine + ".dat".
 An alternative output file can be named as:
    ^blah.dat       send output to "blah.dat" 
 
COMPILATION:

    gcc -lm -lgsl -lgslcblas patchsf.c -o patchSF

This code depends on the installation of Gnu Scientific Library (GSL):

    http://www.gnu.org/software/gsl/

HISTORY:

    Written by Curtis J. Saxton (2011,2012).
    Tidying by Curtis J. Saxton (20160306).
    Updated "README" by Curtis J. Saxton (20160513).

If you apply or adapt the structure function code in academic research, please cite the original method paper as a reference:

    "Long-term X-ray variability of Swift J1644+57"
     Curtis J. Saxton,  Roberto Soria,  Kinwah Wu,  N. Paul M. Kuin,
     2012, Monthly Notices Royal Astronomical Society 422, 1625-1639.
     http://dx.doi.org/10.1111/j.1365-2966.2012.20739.x
     DOI: 10.1111/j.1365-2966.2012.20739.x
     http://arxiv.org/abs/1201.5210
     
