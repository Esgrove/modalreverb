# Modal Filter Reverb

![alt text](https://github.com/Esgrove/modalreverb/blob/master/waterfall.png)

I wrote a seminar paper on modal filter reverberation for the Aalto University Acoustics and Audio Technology Seminar 2016 course, and implemented the reverb algorithm described in the paper [A modal architecture for artificial reverberation with application to room acoustics modeling](http://www.aes.org/e-lib/browse.cfm?elib=17531). In the modal filter reverb audio effect, artificial reverberation is produced by modeling a space (or other reverb producing object, such as analog spring / plate reverb) as a linear combination of resonant filters, where each filter corresponds to one (room) mode of vibration. This way, a digital reverb audio effect can be produced efficiently and accurately from any acoustic space or object. As each mode is independend from all others, the calculation can be easily parallelized and implemented also on a GPU.

### Included in the repo:

 - Seminar paper (Modal Filter Reverberation.pdf)
 - Presentation  (Presentation.pdf)
 - MATLAB scripts for the reverb and images for the paper
 - Audio examples

#### Modal Reverb implementation (modalreverb.m) 

Testing the reverb algorithm with different numbers of filters (modeled modes). Modes are divided randomly inside each octave band, with the number of modes per octave increasing exponentially (2x per each octave). Some randomness is introduced also for the reverberation time as a function of frequency. Parameters are completely synthetic (though following a realistic distribution) so the end result is not the most realistic or best sounding reverb. The intention here was to see how it works and primarily to examine how the number of modes affects the result. *2018 update:* Code cleaned and improved. Now vectorized properly and unnecessary steps removed -> reverb filter loop runs much faster, and the majority of the execution time of the whole script is spent on drawing and exporting the figures.

#### Room frequency response plotting (roomresponse.m)

Plotting a regular 2d frequency response and a fancy waterfall plot from a room acoustic measurement (my living room with Genelec 6010's in this case). The frequency response is plotted directly from measurement data containing the frequencies and magnitudes, whereas the waterwall plot is calculated from the original impulse response. 

#### Room modes (roommodes.m)

Calculate and plot the count and distribution of room modes in an ideal room for reference.