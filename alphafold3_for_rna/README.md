# Has AlphaFold 3 reached its success for RNAs?

This repository contains the code discussed in the Medium article : "Has AlphaFold 3 reached its success for RNAs?".

![](img/best_worst_challenges.gif)


## Installation

To install the libraries needed for the visualisation, you just need to install `plotly` and `sklearn` using:

```bash
pip install plotly scikit-learn
```

## Data

We provide the structures for the five different datasets in `.pdb` format in the `data/pdb` folder. It is composed of the native structures and the predictions from ten different models as well as AlphaFold 3. 

Here is an example of the different predictions for the challenge R1107 of CASP-RNA:

![](img/R1107_preds.gif)

The metrics computed on these predictions are available in the `data/output` folder. 

There is also the `data/plots` folder where are saved the results of the visualisation. 

## Commands

To run the visualisation, you can use the following command:

```bash
python -m src.bar_helper
```

It will output and save the polar chart in `data/plots` folder.


![](img/results_viz.png)


