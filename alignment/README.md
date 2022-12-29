# Alignment Code

This folder contains utility code relevant for physical (optical) alignment as well as the code suite for running spatial calibrations (alignSLMtoCam).

* `./alignment_tools.m` -> code snippets for basic control of various rig components useful during alignment (eg. blanking the SLM, sending holograms to SLM, capture simple image on basler)
* `./getPSF.m` -> runs a z-stack to capture radial and axial PSF of hologram
* `./alignSLMtoCam` -> 3D spatial calibration code


### To Do:
- [ ] refactor getPSF to use new basler and sutter ifaces
- [ ] check new basler for exposure level error


## AlignSLMtoCam
