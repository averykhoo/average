"""
direct python port of logmean.c (revision 191, Modified Fri Jul 13 10:18:28 2012 UTC by sund)
"""
import math
from functools import lru_cache
from typing import List

NMAX = 200000


#
# static SURVO_DATA dat;
# static int xvar;
# static int results_line;
# static int prind;
# static double x[NMAX],logx[NMAX];
# static int n;
# static double lmean;
# static int method; // 1=direct 2=series expansion 3=series expansion (recursive)
#                    // 4=(n-1)th diff.ratio
# static int term_comp;
# static double mean,geom_mean,hmean;
# static double scale_factor;
# static double *d1,*d2;
#
# static int load_data();
# static int comp1();
# static int comp_x(int k);
# static int comp4();
# static int comp5();
# static int comp2();
# static double polm(int n,int m);
# static int comp3();
# static int next_m_distr(int n,int m,int *elem1);
# static int other_means();
# static int print_line(char *x);
#
# /**********************************
# char *specs0[]={ "VARS", "MASK", "IND", "CASES", "SELECT",
#                  "METHOD", "TERM_COMP", "OTHERS", "POWMAX",
#                  "RESULTS", "PRIND",  "!" };
# char **specs=specs0;
# **************************************/
# /****************************
# void main(argc,argv)
# int argc; char *argv[];
# ****************************/
#
# void muste_logmean(char *argv)
#         {
#         int i;
#
# //      if (argc==1) return;
#         s_init(argv);
#
#         if (g<2)
#             {
#             sur_print("\nUsage: LOGMEAN <data>,<var>,L");
#             WAIT; return;
#             }
#         results_line=0;
#         if (g>3)
#             {
#             results_line=edline2(word[3],1,1);
#             if (results_line==0) return;
#             }
#
#     d1=NULL;
#     d2=NULL;
#
#         i=spec_init(r1+r-1); if (i<0) return;
#         i=data_read_open(word[1],&dat);
#         if (i<0) return;
#
#         i=conditions(&dat); if (i<0) return;
#
#         prind=0;
#         i=hae_apu("prind",sbuf); if (i) prind=atoi(sbuf);
#         if ((i=spfind("PRIND"))>=0) prind=atoi(spb[i]);
#
#         xvar=varfind(&dat,word[2]); if (xvar<0) return;
#
#         method=3;
#         i=spfind("METHOD");
#         if (i>=0) method=atoi(spb[i]);
#
#         term_comp=0;
#         i=spfind("TERM_COMP"); if (i>=0) term_comp=atoi(spb[i]);
#
#         i=load_data(); if (i<0) return;
#         if (n<2L) { sur_print("\nAt least 2 observations needed!");
#                     WAIT; return;
#                   }
#
#         if (method>3)
#             {
#             if (method==5)
# //              for (i=0; i<n; ++i)
# //              logx[i]=sqrt(x[i]);
#             comp5();
#
#             else comp4();
#             }
#         else if (method==3) comp3();
#         else if (method==2) comp2();
#         else if (method==1) comp1();
#         else comp_x(-method);
#
#         i=output_open(eout); if (i<0) return;
#         sprintf(sbuf,"Data: %s Variable: %s  N=%d",word[1],word[2],n);
#         print_line(sbuf);
#         sprintf(sbuf,"Logarithmic mean: %16.16g",lmean/scale_factor);
#         print_line(sbuf);
#         i=spfind("OTHERS");
#         if (i>=0 && atoi(spb[i])>0)
#             {
#             other_means();
#             sprintf(sbuf,"Arithmetic mean:  %16.16g",mean/scale_factor);
#             print_line(sbuf);
#             sprintf(sbuf,"Geometric mean:   %16.16g",geom_mean/scale_factor);
#             print_line(sbuf);
#             sprintf(sbuf,"Harmonic mean:    %16.16g",hmean/scale_factor);
#             print_line(sbuf);
#             }
#         output_close(eout);
#         s_end(argv);
#         }
#

def muste_logmean(argv):
    i = 0

    s_init(argv)

    if g < 2:
        sur_print("\nUsage: LOGMEAN <data>,<var>,L")
        WAIT
        return

    results_line = 0
    if g > 3:
        results_line = edline2(word[3], 1, 1)
        if results_line == 0:
            return

    d1 = None
    d2 = None

    i = spec_init(r1 + r - 1)
    if i < 0:
        return

    i = data_read_open(word[1], & dat)
    if i < 0:
        return

    i = conditions( & dat)
    if i < 0:
        return

    prind = 0
    i = hae_apu("prind", sbuf)
    if i:
        prind = atoi(sbuf)
    if ((i = spfind("PRIND")) >= 0):
        prind = atoi(spb[i])

    xvar = varfind( & dat, word[2])
    if xvar < 0:
        return

    method = 3
    i = spfind("METHOD")
    if i >= 0:
        method = atoi(spb[i])

    term_comp = 0
    i = spfind("TERM_COMP")
    if i >= 0:
        term_comp = atoi(spb[i])

    i = load_data()
    if i < 0:
        return
    if n < 2:
        sur_print("\nAt least 2 observations needed!")
        WAIT
        return

    if method > 3:
        if method == 5:
            comp5()
        else:
            comp4()
    elif method == 3:
        comp3()
    elif method == 2:
        comp2()
    elif method == 1:
        comp1()
    else:
        comp_x(-method)

    i = output_open(eout)
    if i < 0:
        return

    i = spfind("OTHERS")
    if i >= 0 and atoi(spb[i]) > 0:
        other_means()

    output_close(eout)
    s_end(argv)


def load_data(data, method, term_comp):
    a: float
    _min: float = 1e308
    _max: float = 1e-308
    x = []
    log_x = []
    for a in data:
        if a <= 0:
            raise ValueError(f'Only positive values allowed! Got {a}')
        if a < _min:
            _min = a
        if a > _max:
            _max = a
        x.append(a)
        log_x.append(math.log(a))
    n = len(x)
    if _min != _max:
        scale_factor = 2 / (_max - _min)
    else:
        scale_factor = _max

    if method <= 0 or term_comp:
        scale_factor = 1

    for l in range(len(x)):
        x[l] *= scale_factor
        log_x[l] = math.log(x[l])

    return n, x, log_x, scale_factor


# static int comp1()
#     {
#     int i,j;
#     double fact;
#     double term;
#     double lx;
#
#     fact=1.0; lmean=0;
#     for (i=0; i<n; ++i)
#         {
#         if (i) fact*=(double)i;
#         term=x[i]; lx=logx[i];
#         for (j=0; j<n; ++j)
#             {
#             if (i!=j) term/=lx-logx[j];
#             }
#         lmean+=term;
#         }
#     lmean*=fact;
#     return(1);
#     }
def comp1(n, x, log_x):
    fact = 1.0
    lmean = 0
    for i in range(n):
        if i:
            fact *= i
        term = x[i]
        lx = log_x[i]
        for j in range(n):
            if i != j:
                term /= lx - log_x[j]
        lmean += term
    lmean *= fact
    return lmean


# static int comp_x(int k)
#     {
#     int i,j;
#     double term;
#     double lx;
#
#     lmean=0;
#     for (i=0; i<n; ++i)
#         {
#         term=pow(x[i],(double)k);
#         lx=x[i];
#         for (j=0; j<n; ++j)
#             {
#             if (i!=j) term/=lx-x[j];
#             }
#         lmean+=term;
#         }
#     return(1);
#     }
def comp_x(n, x, log_x, k):
    lmean = 0
    for i in range(n):
        term = math.pow(x[i], k)
        lx = x[i]
        for j in range(n):
            if i != j:
                term /= lx - x[j]
        lmean += term
    return lmean


# static int comp4()
#     {
#     int i,j;
#     double fact;
#
#
#
#     d1=(double *)muste_malloc(n*sizeof(double));
#     d2=(double *)muste_malloc(n*sizeof(double));
#
#     fact=1.0; lmean=0;
#     for (i=1; i<n; ++i)
#         fact*=(double)i;
#
#     for (i=0; i<n; ++i) d1[i]=x[i];
#
# // for (i=0; i<n; ++i) Rprintf("\n%d %g %g",i,d1[i],logx[i]);
#
#     for (j=0; j<n-1; ++j)
#         {
#         for (i=0; i<n-j-1; ++i)
#             {
# // Rprintf("\nj=%d i=%d %g %g %g %g",j,i,d1[i+1],d1[i],logx[i+j+1],logx[i]);
# // getch();
#             d2[i]=(d1[i+1]-d1[i])/(logx[i+j+1]-logx[i]);
#             }
#         for (i=0; i<n-j-1; ++i) d1[i]=d2[i];
#         }
#     lmean=fact*d1[0];
#     return(1);
#     }
def comp4():
    i = 0
    j = 0
    fact = 0

    d1 = (double *)(muste_malloc(n * sizeof(double)))
    d2 = (double *)(muste_malloc(n * sizeof(double)))

    fact = 1.0
    lmean = 0
    for i in range(1, n):
        fact *= (double)(i)

    for i in range(0, n):
        d1[i] = x[i]

    for j in range(0, n - 1):
        for i in range(0, n - j - 1):
            d2[i] = (d1[i + 1] - d1[i]) / (logx[i + j + 1] - logx[i])
        for i in range(0, n - j - 1):
            d1[i] = d2[i]
    lmean = fact * d1[0]
    return 1


# static int comp5()
#     {
#     int i,j;
#     double fact;
#
#
#
#     d1=(double *)muste_malloc(n*sizeof(double));
#     d2=(double *)muste_malloc(n*sizeof(double));
#
#     fact=1.0; lmean=0;
# //  for (i=1; i<n; ++i)
# //      fact*=(double)i;
#
#     for (i=0; i<n; ++i) d1[i]=x[i];
#
# // for (i=0; i<n; ++i) Rprintf("\n%d %g %g",i,d1[i],logx[i]);
#
#     for (j=0; j<n-1; ++j)
#         {
#         for (i=0; i<n-j-1; ++i)
#             {
# // Rprintf("\nj=%d i=%d %g %g %g %g",j,i,d1[i+1],d1[i],logx[i+j+1],logx[i]);
# // getch();
#
#             d2[i]=fact*(d1[i+1]-d1[i])/(logx[i+j+1]-logx[i]);
#             }
#         for (i=0; i<n-j-1; ++i) d1[i]=d2[i];
#         ++fact;
#         }
#     lmean=d1[0];
#     return(1);
#     }
#
def comp5():
    i = 0
    j = 0
    fact = 0

    d1 = (double *)(muste_malloc(n * sizeof(double)))
    d2 = (double *)(muste_malloc(n * sizeof(double)))

    fact = 1.0
    lmean = 0

    for i in range(0, n):
        d1[i] = x[i]

    for j in range(0, n - 1):
        for i in range(0, n - j - 1):
            d2[i] = fact * (d1[i + 1] - d1[i]) / (logx[i + j + 1] - logx[i])
        for i in range(0, n - j - 1):
            d1[i] = d2[i]
        fact += 1

    lmean = d1[0]
    return 1


# #define MMAX 200000
# #define POWMAX 60
#
# static int pot[MMAX];
# static double powlog[MMAX][POWMAX];
# static int powmax=POWMAX;
# static double pm[MMAX][POWMAX];
#


# static int comp2()
#     {
#     int i;
#     int m;
#     int ncomb;
#     double term1,term2;
#     double fact;
#     double lmean1=0.0;
#
#     i=spfind("POWMAX");
#     if (i>=0) powmax=atoi(spb[i]);
#     if (powmax>POWMAX) powmax=POWMAX;
#
#     for (i=0; i<n; ++i) powlog[i][0]=logx[i];
#     lmean=1.0;
#     m=1; fact=1.0;
#     while (1)
#         {
# // Rprintf("\nm=%d",m); getch();
#         for (i=0; i<n; ++i) pot[i]=0;
#         pot[n-1]=m; ncomb=1L;
#         term2=0.0;
#         while (1)
#             {
# //          sur_print("\n");
# //          for (i=0; i<n; ++i) { sprintf(sbuf,"%d ",pot[i]); sur_print(sbuf); }
#             term1=1.0;
#             for (i=0; i<n; ++i)
#                 {
#                 if (pot[i]!=0) term1*=powlog[i][pot[i]-1];
#                 }
#             term2+=term1;
# // sprintf(sbuf,"\nterm1=%g term2=%g",term1,term2); sur_print(sbuf); getch();
#             i=next_m_distr(m,n,pot);
#             if (i<0) break;
#             ++ncomb;
#             }
#         fact*=(double)m;
# //      Rprintf("\nncomb=%ld fact=%g",ncomb,fact); getch();
#         lmean+=term2/(double)ncomb/fact;
# sprintf(sbuf,"\n%d: lmean=%16.16g",m,lmean); sur_print(sbuf);
#         if (fabs(lmean-lmean1)==0.0 || m>=powmax) break;
#         lmean1=lmean;
#         ++m;
#
#         for (i=0; i<n; ++i) powlog[i][m-1]=logx[i]*powlog[i][m-2];
#         }
#     return(1);
#     }
#
POWMAX = 100  # Maximum value of the powmax variable

# Initialize the powlog array
powlog = [[0 for _ in range(POWMAX)] for _ in range(n)]


def comp2():
    # Read the value of the powmax variable
    powmax = POWMAX
    i = spfind("POWMAX")
    if i >= 0:
        powmax = atoi(spb[i])
    if powmax > POWMAX:
        powmax = POWMAX

    # Initialize the powlog and pot arrays
    for i in range(n):
        powlog[i][0] = logx[i]
    lmean = 1.0
    m = 1
    fact = 1.0
    while True:
        for i in range(n):
            pot[i] = 0
        pot[n - 1] = m
        ncomb = 1
        term2 = 0.0
        while True:
            term1 = 1.0
            for i in range(n):
                if pot[i] != 0:
                    term1 *= powlog[i][pot[i] - 1]
            term2 += term1
            i = next_m_distr(m, n, pot)
            if i < 0:
                break
            ncomb += 1
        fact *= m
        lmean += term2 / ncomb / fact
        if fabs(lmean - lmean1) == 0.0 or m >= powmax:
            break
        lmean1 = lmean
        m += 1

        for i in range(n):
            powlog[i][m - 1] = logx[i] * powlog[i][m - 2]

    return 1


# static double polm(int n,int m)
#     {
#     int i;
#     double s;
#
#     s=pm[n-1][m-1];
#     if (s!=1e308) return(s);
#
#     if (m==1)
#         {
#         s=0.0;
#         for (i=0; i<n; ++i) s+=logx[i];
#         pm[n-1][m-1]=s;
#         return(s);
#         }
#     if (n==1)
#         {
#         s=pm[n-1][m-1]=powlog[0][m-1];
#         return(s);
#         }
#     s=powlog[n-1][m-1];
#     for (i=1; i<m; ++i) s+=powlog[n-1][m-i-1]*polm(n-1,i);
#     s+=polm(n-1,m);
#     pm[n-1][m-1]=s;
#     return(s);
#     }

# output by chatgpt

# Initialize the pm and powlog arrays
pm = [[math.inf for _ in range(m)] for _ in range(n)]
powlog = [[0 for _ in range(m)] for _ in range(n)]


def polm(n: int, m: int) -> float:
    # If pm[n-1][m-1] has been calculated before, return that value
    if pm[n - 1][m - 1] != math.inf:
        return pm[n - 1][m - 1]

    # If m equals 1, sum the log of each x value in the range from 0 to n
    if m == 1:
        s = sum(logx[:n])
        pm[n - 1][m - 1] = s
        return s

    # If n equals 1, get the m-1th power of the log of the 0th element of powlog
    if n == 1:
        pm[n - 1][m - 1] = powlog[0][m - 1]
        return pm[n - 1][m - 1]

    # Calculate the mth power of the log of the nth element of powlog
    # and add it to the sum of the (m-1)th power of the log of the (n-1)th element of powlog
    # multiplied by the result of calling polm with n-1 and i, for each i in the range from 1 to m
    # Then add the result of calling polm with n-1 and m
    s = powlog[n - 1][m - 1]
    for i in range(1, m):
        s += powlog[n - 1][m - i - 1] * polm(n - 1, i)
    s += polm(n - 1, m)

    pm[n - 1][m - 1] = s
    return s


#
# given the following code:
#
# ```
# for (i=0; i<n; ++i) powlog[i][0]=logx[i];
# ...
# for (i=0; i<n; ++i)
#     {
#     if (pot[i]!=0) term1*=powlog[i][pot[i]-1];
#     }
# ...
# for (i=0; i<n; ++i) powlog[i][m-1]=logx[i]*powlog[i][m-2];
# ```
# describe powlog

# powlog is a two-dimensional array where each element powlog[i][j] represents the jth power of the natural logarithm
# of the ith element of x. It is used to store the results of these calculations so they can be used later without
# having to recompute them.
#
# In the first line of the code snippet, powlog[i][0] is set to logx[i], which is the natural logarithm of the
# ith element of x.
#
# In the second line of the code snippet, powlog[i][pot[i]-1] is multiplied with term1 for each i in the range
# from 0 to n. pot[i] is assumed to be the power that the ith element of x is raised to, so powlog[i][pot[i]-1]
# is the value of the (pot[i]-1)th power of the natural logarithm of the ith element of x.
#
# In the third line of the code snippet, powlog[i][m-1] is set to the product of logx[i] and powlog[i][m-2].
# This calculates the m-1th power of the natural logarithm of the ith element of x by multiplying the (m-2)th power
# by the natural logarithm of x[i].
#
#
#
#
# averykhoo@gmail.com
# logx is equivalent to [math.log(x) for x in xs]. describe polm

# polm() is a recursive function that calculates the mth power of the natural logarithm of the sum of the elements
# in a list x. It uses a two-dimensional array pm as a cache to store previously calculated results, so that they
# can be re-used instead of being recomputed.
#
# The function takes two arguments: n, the length of the list x, and m, the power to which the natural logarithm
# of the sum of the elements in x should be raised. If m equals 1, the function returns the sum of the natural
# logarithm of each element in x. If n equals 1, the function returns the m-1th power of the natural logarithm of
# the first element in x. Otherwise, the function calculates the mth power of the natural logarithm of the sum of
# the elements in x by summing the mth power of the natural logarithm of each element in x, the (m-1)th power of
# the natural logarithm of each element in x multiplied by the result of calling polm() with n-1 and i for each i
# in the range from 1 to m, and the result of calling polm() with n-1 and m. The result is then cached in pm and
# returned.

@lru_cache
def polm_cached(n, m):
    # If m equals 1, sum the log of each x value in the range from 0 to n
    if m == 1:
        return sum(logx[:n])

    # If n equals 1, get the m-1th power of the log of the 0th element of powlog
    if n == 1:
        return powlog[0][m - 1]

    # Calculate the mth power of the log of the nth element of powlog
    # and add it to the sum of the (m-1)th power of the log of the (n-1)th element of powlog
    # multiplied by the result of calling polm with n-1 and i, for each i in the range from 1 to m
    # Then add the result of calling polm with n-1 and m
    s = powlog[n - 1][m - 1]
    for i in range(1, m):
        s += powlog[n - 1][m - i - 1] * polm(n - 1, i)
    s += polm(n - 1, m)

    return s


from functools import lru_cache


@lru_cache
def polm2(xs, m):
    # If m equals 1, sum the log of each element in xs
    if m == 1:
        return sum(math.log(x) for x in xs)

    # If the length of xs equals 1, get the m-1th power of the log of the 0th element of xs
    if len(xs) == 1:
        return math.pow(math.log(xs[0]), m - 1)

    # Calculate the mth power of the log of the sum of the elements in xs
    # and add it to the sum of the (m-1)th power of the log of each element in xs
    # multiplied by the result of calling polm with the sublist of xs without the first element
    # and i, for each i in the range from 1 to m
    # Then add the result of calling polm with the sublist of xs without the first element and m
    s = math.pow(math.log(sum(xs)), m)
    for i in range(1, m):
        s += math.pow(math.log(xs[i]), m - 1) * polm(xs[1:], i)
    s += polm(xs[1:], m)

    return s


# static int comp3()
#     {
#     int i;
#     int m;
#     double ncomb;
#     double term2;
#     double fact;
#     double lmean1=0.0;
#     double term3,aterm,gterm;
#     double sum=0;
#
#
#     i=spfind("POWMAX");
#     if (i>=0) powmax=atoi(spb[i]);
#     if (powmax>POWMAX) powmax=POWMAX;
#
#     if (term_comp)
#         {
#         sum=0.0;
#         for (i=0; i<n; ++i) sum+=logx[i];
#         sum/=(double)n;
#         }
#
#     for (i=0; i<n; ++i) for(m=0; m<powmax; ++m) pm[i][m]=1e308;
#
#     for (i=0; i<n; ++i)
#         {
#         powlog[i][0]=logx[i];
#         for (m=2; m<=powmax; ++m) powlog[i][m-1]=powlog[i][m-2]*logx[i];
#         }
#
#     lmean=1.0;
#     m=1; fact=1.0; ncomb=n;
#     while (1)
#         {
# // Rprintf("\nm=%d",m); getch();
#         term2=polm(n,m);
# // Rprintf("\nterm2=%g",term2); getch();
#         fact*=(double)m;
# //      Rprintf("\nncomb=%g fact=%g",ncomb,fact); getch();
#         if (!term_comp)
#             lmean+=term2/(double)ncomb/fact;
#         else
#             {
#             term3=term2/(double)ncomb;
#             lmean+=term2/(double)ncomb/fact;
#             aterm=0.0;
#             for (i=0; i<n; ++i) aterm+=pow(logx[i],(double)m);
#             aterm/=(double)n;
#             gterm=pow(sum,(double)m);
#             if (term3<gterm || aterm<term3)
#    {
# sprintf(sbuf,"\nm=%d gterm=%g term3=%g aterm=%g",m,gterm,term3,aterm);
# sur_print(sbuf); sur_getch();
#    }
#
#             }
# // sprintf(sbuf,"\n%d: lmean=%16.16g",m,lmean); sur_print(sbuf);
#         if (fabs(lmean-lmean1)==0.0 || m>=powmax) break;
#         lmean1=lmean;
#
#         ++m;
#         ncomb*=(double)(n+m-1)/(double)m;
#         }
#     return(1);
#     }
#

POWMAX = 100  # Maximum value of the powmax variable

# Initialize the pm and powlog arrays
pm = [[inf for _ in range(powmax)] for _ in range(n)]
powlog = [[0 for _ in range(powmax)] for _ in range(n)]


def comp3():
    # Read the value of the powmax variable
    powmax = POWMAX
    i = spfind("POWMAX")
    if i >= 0:
        powmax = atoi(spb[i])
    if powmax > POWMAX:
        powmax = POWMAX

    # Initialize the sum variable
    sum = 0.0
    if term_comp:
        for i in range(n):
            sum += logx[i]
        sum /= n

    # Initialize the pm and powlog arrays
    for i in range(n):
        for m in range(powmax):
            pm[i][m] = inf
        powlog[i][0] = logx[i]
        for m in range(2, powmax + 1):
            powlog[i][m - 1] = powlog[i][m - 2] * logx[i]

    lmean = 1.0
    m = 1
    fact = 1.0
    ncomb = n
    while True:
        term2 = polm(n, m)
        fact *= m
        if not term_comp:
            lmean += term2 / ncomb / fact
        else:
            term3 = term2 / ncomb
            lmean += term2 / ncomb / fact
            aterm = 0.0
            for i in range(n):
                aterm += pow(logx[i], m)
            aterm /= n
            gterm = pow(sum, m)
            if term3 < gterm or aterm < term3:
                # Handle this case
                pass
        if fabs(lmean - lmean1) == 0.0 or m >= powmax:
            break
        lmean1 = lmean
        m += 1
        ncomb *= (n + m - 1) / m

    return 1


# static int next_m_distr(int n,int m,int *elem1) // n ja m vaihtaneet paikkojaan!
#         {
#         int i,k;
#
#         i=m-1;
#         while (elem1[i]==0) --i;
#         if (i==0) return(-1);
#         ++elem1[i-1];
#         elem1[m-1]=elem1[i]-1;
#         for (k=i; k<m-1; ++k) elem1[k]=0;
#         return(1);
#         }
#

def next_m_distr(n: int, m: int, elem1: List[int]) -> int:
    """
    This function generates the next element of a distribution with m elements in n classes.
    The n and m parameters are in the opposite order compared to the original C code!
    """
    i = m - 1
    while elem1[i] == 0:
        i -= 1
    if i == 0:
        return -1
    elem1[i - 1] += 1
    elem1[m - 1] = elem1[i] - 1
    for k in range(i, m - 1):
        elem1[k] = 0
    return 1


# static int other_means()
#     {
#     long l;
#
#     mean=geom_mean=hmean=0.0;
#     for (l=0L; l<n; ++l)
#         {
#         mean+=x[l]; geom_mean+=logx[l];
#         hmean+=1.0/x[l];
#         }
#     mean/=(double)n; geom_mean=exp(geom_mean/(double)n);
#     hmean=(double)n/hmean;
#     return(1);
#     }
#
#
def other_means(xs):
    n = len(xs)
    logx = [math.log(x) for x in xs]
    mean = geom_mean = hmean = 0.0
    for l in range(n):
        mean += xs[l]
        geom_mean += logx[l]
        hmean += 1.0 / xs[l]
    mean /= n
    geom_mean = math.exp(geom_mean / n)
    hmean = n / hmean
    return mean, geom_mean, hmean

# static int print_line(char *x)
#         {
#         output_line(x,eout,results_line);
#         if (results_line) ++results_line;
#         if (results_line>r2) results_line=0;
#         return(1);
#         }
