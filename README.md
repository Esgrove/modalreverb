# Modal Filter Reverb

![alt text](https://github.com/Esgrove/modalreverb/blob/master/waterfall.png)

I wrote a seminar paper on modal filter reverberation for the Aalto University Acoustics and Audio Technology Seminar 2016 course, and implemented the reverb algorithm described in the paper [A modal architecture for artificial reverberation with application to room acoustics modeling](http://www.aes.org/e-lib/browse.cfm?elib=17531). In a modal filter reverb, artificial reverberation is produced by modeling a space (or other reverb producing object, such as analog spring reverbs) as a linear combination of resonant filters, where each filter corresponds to one (room) mode. This way, a digital reverb audio effect can be produced efficiently and accurately from any acoustic space or object.

### Included in the repo:

 - Seminar paper (Modal Filter Reverberation.pdf)
 - Presentation  (Presentation.pdf)
 - MATLAB scripts for the reverb and various images used in the paper:


#### Modal Reverb implementation (modalreverb.m) 

Testing the reverb algorithm with different numbers of filters (modeled modes). Modes are divided randomly inside each octave band, with the number of modes per octave increasing exponentially. Some randomness is introduced also for the reverberation time as a function of frequency. Parameters are completely synthetic so the end result is not the most realistic or best sounding. The intention here was to see how it works and mainly see how the number of modes affects the result. The reverb code is far from optimal and is quite slow for long audio samples (it is programmed in a very straigthforward manner, i.e. using for-loops instead of using MATLAB vectorization).

#### Room frequency response plotting (roomresponse.m)

Plotting a regular 2d frequency response and a fancy waterfall plot from a room acoustic measurement (my living room in this case). The frequency response is plotted directly from measurement data containing the frequencies and magnitudes, whereas the waterwall plot is calculated from the original impulse response. 

#### Room modes (roommodes.m)

Calculate and plot the count and distribution of room modes in a room for reference.