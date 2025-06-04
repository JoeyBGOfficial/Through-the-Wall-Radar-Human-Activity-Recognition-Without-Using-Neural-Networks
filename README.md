# Through-the-Wall-Radar-Human-Activity-Recognition-Without-Using-Neural-Networks
## I. Introduction ##

**Write Sth. Upfront:** 

This paper is dedicated to the memory of my grandma.

This is probably the only paper in my life that was written in a hospital the entire time. This is also at the same time probably the only work in my life that is publicly available on preprint platform arXiv, not submitted to a formal journal, but open-sourced. It's a very unique idea. Think of it as my last time of crazy burn for the field I've been fighting for nearly $5$ years during my PhD carrier.

We have been stuck in using neural network models to achieve radar target recognition for so many years. I just want to be back to the 90s and 00s, when people could also achieve complex tasks with a certain level of intelligence with perfect physical interpretability using traditional and recent signal processing techniques. It was this vision that inspired me and this work was born.

I want to thank my grandpa, mother, father, aunt, for the impeccable care you gave to my grandmother on her deathbed and for making me feel the most precious affection on earth. I want to thank my friends who cared for my family during this time, and my mentors for nurturing and trusting me with the ability to accomplish this work. Additionally, thanks to my love Miss Xu, it is our commitment that makes this work possible.

I truly hope everyone can find something in my journey. It is the power of love that penetrates all difficulties.

![Introduction](https://github.com/user-attachments/assets/1fdef49f-98f2-4b03-ac26-fef91f58b39c)

Fig. 1. Current works in this field take neural network-based methods as the research hotspot. This work returns to rethink the value of traditional mindsets.

**Basic Information:** This repository is the open source code for my latest work: "Through-the-Wall Radar Human Activity Recognition WITHOUT Using Neural Networks", submitted to arXiv.

**Email:** JoeyBG@126.com;

**Abstract:** After a few years of research in the field of through-the-wall radar (TWR) human activity recognition (HAR), I found that we seem to be stuck in the mindset of training on radar image data through neural network models. The earliest related works in this field based on template matching did not require a training process, and I believe they have never died. Because these methods possess a strong physical interpretability and are closer to the basis of theoretical signal processing research. In this paper, I would like to try to return to the original path by attempting to eschew neural networks to achieve the TWR HAR task and challenge to achieve intelligent recognition as neural network models. In detail, the range-time map and Doppler-time map of TWR are first generated. Then, the initial regions of the human target foreground and noise background on the maps are determined using corner detection method, and the micro-Doppler signature is segmented using the multiphase active contour model. The micro-Doppler segmentation feature is discretized into a two-dimensional point cloud. Finally, the topological similarity between the resulting point cloud and the point clouds of the template data is calculated using Mapper algorithm to obtain the recognition results. The effectiveness of the proposed method is demonstrated by numerical simulated and measured experiments.

**Corresponding Papers:**

[1]
