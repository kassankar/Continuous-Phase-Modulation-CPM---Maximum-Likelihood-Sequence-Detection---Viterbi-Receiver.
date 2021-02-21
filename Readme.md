# Continuous Phase Modulation (CPM) - Maximum Likelihood Sequence Detection - Viterbi Receiver.

A Matlab code for the "Maximum Likelihood Sequence Detection" for any CPM modulation (GMSK-RECT-RC....) using the Viterbi Algorithm. This work is based on the book: __' Digital Communication (Proakis)'__, and thesis: __Comparison of Noncoherent detectors for SOQPSK and GMSK in Phase Noise Channels__.


# How to cite this Work ?



# How to run the code ?
1. Make sure that you have a compatible version of Matlab (this code was tested on Matlab 2018b)
2. Download (clone) the files from the repository.
3. Open the file called _MLSD_Viterbi_CPM.m_
4. Select the section called __Pulse shape & Variable ini__
5. Select the type of pulses by changing the variable `pulse` number
    * `1` is for Lorentzian
    * `2` is for GMSK
    * `3` is for Raised Cosine
    * `4` is for Rectangular
	* We can define a new pulse if the wanted pulse is not presented. Go to the file called __MainFunctions.m__, in the first function __CREATECPMPULSE__ and add the pulse, using the same method used to define the other pulses.
6. Change the pulse length by changing the variable `pulse_length`.
7. Select the oversampling by changing the variable `os`.
8. Select M-ary by changing the `M_ary` (e.g., M=2 for Binary).
9. Select the minimum Euclidean distance of the modulation by changing the variable `dmin`.
10. Select after how many observation symbols the Viterbi receiver can make a decision on the detected bit by changing the variable `decision_delay` (usually 50 is high enough). 
## Note:
For points number `9` and `10`. We can obtain `dmin` and the number of observation symbols from the code we provided in this link: https://github.com/kassankar/ContinuousPhaseModulation-CPM---EuclideanDistance.


# Example

In this example we present the BER of __GMSK__ CPM modulation with `BT=0.3`. 
````
pulse            = 2;
pulse_length     = 3;     
os               = 2^2;    
Ts               = 1/os;   
M_ary            = 2^1;    
modulation_index = 0.5;    
dmin             = 1.78;   % GMSK BT=0.3
decision_delay   = 50;
````
Results Plot:
![](BER_GMSK.png?raw=true)
## Known bug
For `pulse_length = 1` (Full response), with `M_ary >2`, we may obtain an error for some modulation indices.
## Warning

## Contributing ?


## License
© 2020-2021 Karim Kassan

