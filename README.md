# A Cone-Beam Computed Tomography (CBCT) Analysis Pipeline

## 🌐 About Me
Hi! I'm Nathan Cao, a current junior at Rice University. I am an electrical and computer enginering major, and outside of academics I pursue medical imaging research at MD Anderson Cancer Center under the guidance of Dr. Ke Li. Of course, research isn't my only extracurricular endeavor. I am a part-time triathlete (hoping to compete at the upcoming Ironman 70.3 Waco), ex-competitive badminton player, and a coffee addict (love a good flat white). Post-graduation I hope to pursue a Ph.D. in medical imaging and one day develop low-radiation imaging modalities suitable for underserved regions. 

---

## 🧠 Problem Statement
CT is one of the most widely used diagnostic and intervention imaging tools. From trauma to brain surgeons, CT is the modality that most rely on. Why? It's simple. It's inexpensive (making it a crowd pleaser among the general population), and fast (so anyone experiencing trauma can get scanned in a timely manner). The project first started when I commenced my research journey with Dr. Li. I had partaken in imaging research before, but because my primary duties involved building simulations, I never got a taste of learning about the field from a classical approach. Dr. Li very generously devoted his time every week to guiding me up the steep learning curve of everything CT, and this project is the result. I hope whoever stumbles across this (very heavily commented) repository uses it however they wish, whether to support their own project ambitions or as a guide. 

---

## 📊 Forward and Back Projection
Like I said previously, this project is honestly the culmination of my educational journey with CT, so several scripts are iterations of one another
- cbct_head_phantom_analysis1: my first iteration of attempting forward projection. The code works just fine but is quite slow
- cbct_head_phantom_analysis: my second iteration of attempting forward projection. The code is heavily inspired by Dr. Li's own reconstruction methods. This contains the classic Siddon forward projection, classic FDK reconstruction, C++ accelerated Siddon (code is in siddon_kernel.cpp), and CUDA accelerated FDK reconstruction (currently still under development)

---

## ⚙️ Acceleration
- siddon_kernel.cpp: this cpp kernel accelerates the classic Siddon forward projection which takes around give or take a few hours on my device to just 10 minutes. For those wondering why, cpp runs closer to the CPU than MATLAB so the execution time is shortened by a lot.
- fdk_kernel.cu: this is a CUDA kernel meant to accelerate FDK reconstruction. The projected reconstruction time should be around 5-10 seconds.

## ⬇️ Downloads
- mingw-84 so that Matlab can use Mex to call C++ code
- Visual Studio Build tools so that Matlab can call CUDA
- CUDA toolkit (whatever version is compatible with the specific version of Matlab)
- Visual Studio (so that you can use the full functionality of CUDA)

