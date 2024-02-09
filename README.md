<br/>
<h1 align="center">synMinE</h1>
<p align="center"><i>Code for the paper "Machine Learning-Aided Design and Screening of an Emergent Protein Function in Synthetic Cells".</i></p>
<br/>

<br/>
This repository contains all data and code for the paper "Machine Learning-Aided Design and Screening of an Emergent Protein Function in Synthetic Cells" by Shunshi Kohyama, Béla P. Frohn, Leon Babl and Petra Schwille, published in Nature Communications (doi).
<br/>

<br/>
<p align="center">
    <img width="70%" src="https://github.com/BelaFrohn/syninE/Data/Fig1.png" alt="synMine_i3_Concept">
</p>
<br/>

## Repository Structure
The repository contains five directories: 
- [Data](https://github.com/BelaFrohn/syninE/Data): Containing all sequences, predicted structures, scores and the model weights. 
- [insilicoScoring](https://github.com/BelaFrohn/syninE/insilicoScoring): The pipeline used to score the generated variants <i>in silico</i>. 
- [PostHocAnalysis](https://github.com/BelaFrohn/syninE/PostHocAnalysis): The jupyter notebook and data used to analyse the <i>in silico</i> scores after unblinding the <i>in vitro</i> results. 
- [ProteinGeneration](https://github.com/BelaFrohn/syninE/ProteinGeneration): The code used to geenrate novel MinE-like sequences. 
- [Analysis](https://github.com/BelaFrohn/syninE/Analysis): The code used to analyse <i>in vitro</i> and <i>in vivo</i> data. 


## 📘&nbsp; License
TODO

## ✏️&nbsp; Citation
If you use this code or our pretrained models for your publication, please cite the original paper:
```
@ARTICLE
{TODO,
author={Kohyama, Shunshi and Frohn, Béla P. and Babl, Leon and Schwille, Petra},
journal={Nature Communications},
title={ProtTrans: Towards Cracking the Language of Lifes Code Through Self-Supervised Deep Learning and High Performance Computing},
year={2024},
volume={},
number={},
pages={},
doi={todo}}
```