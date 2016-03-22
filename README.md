# patchSF
Given piecewise patchy time series data, calculate the Kolmogorov structure functions at given order, and respective uncertainties.  The code "patchsf" calculates the Kolmogorov structure function (at given order) from discrete time series data fron a 4-column input file.  Write a 4-column output file containing timescales, SF denominator, SF numerator, uncertainty.

The first command-line argument is the input data file.

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
    
HISTORY:

    Written by Curtis J. Saxton (2011,2012).
    Tidying by Curtis J. Saxton (20160306).

If you apply or adapt the structure function code in academic research, please cite the original method paper as a reference:

    "Long-term X-ray variability of Swift J1644+57"
     Curtis J. Saxton,  Roberto Soria,  Kinwah Wu,  N. Paul M. Kuin,
     2012, Monthly Notices Royal Astronomical Society 422, 1625-1639.
     <a href="http://dx.doi.org/10.1111/j.1365-2966.2012.20739.x">http://dx.doi.org/10.1111/j.1365-2966.2012.20739.x</a>
     http://dx.doi.org/10.1111/j.1365-2966.2012.20739.x
     DOI: 10.1111/j.1365-2966.2012.20739.x
     http://arxiv.org/abs/1201.5210
     
