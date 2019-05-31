This is a model of the ELL circuitry in weakly electric fish
(A. leptorynchus) and how it cancels redundant input using the
indirect feedback pathway (i.e. cerebellar-like granule cells) and a
burst-dependent learning rule as described in this paper:

Bol K, Marsat G, Harvey-Girard E, Longtin A, Maler L (2011)
Frequency-Tuned Cerebellar Channels and Burst-Induced LTD Lead to the
Cancellation of Redundant Sensory Inputs. J Neurosci 31:11028-38

The ELL model is presented in 2 formats here:

In a MATLAB format (CLSglobal.m) with an associated MEX file
(LIFDAPmatlab.c)

And completely in C format (CLSglobal.c, LIFDAPC.c and
BaysDurhamrand.c)

THe MATLAB format was constructed first but the MEX file type (which
requires installing within MATLAB a C compiler and compiling
LIFDAPmatlab.c within MATLAB) can be difficult to implement on some
computers,

To compile LIFDAPmatlab.c within MATLAB IF YOU HAVE A COMPILER IN
MATLAB INSTALLED, type: mex LIFDAPmatlab.c Then run CLSglobal.m

For the C code, one will need to compile the code using a C compiler.

The inputs to both are the frequency of the redundant stimulus and
whether the feedback pathway is active or not (1 or 0, respectively).

The output of both are the firing times as a function of the phase of
the input stimulus (i.e. a post-stimulus time histogram), the
interspike interval histogram up to 200ms, the average firing rate,
the average 2-spike burst rate, the average 4-spike burst rate, and
the average weight value.

Also included (although one will have to uncomment some code) are 2
different ways to model reduced granule cell firing activity after
10Hz - either by decreasing the learning rate, eta, or the feedback
strength, Lambda.


E-mail me at kieran_bol@hotmail.com, or Andre Longtin at
alongtin@uottawa.ca if you have trouble implementing the code.
