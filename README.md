{\rtf1\ansi\ansicpg1252\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww17080\viewh11740\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Code for: \'93Ecological theory helps reveal human decision rules from empirical fishery data\'94\
\
This respiratory contains the simulation and analysis code for the following manuscript:\
> Type the article and authors here\
\
## File structure\
The whole analysis is divided into three parts: abundance distribution prediction, individual-based model, and model performance evaluation, each in its own folder. \
\
The contents in each folder include:\
- `01_abundance_predict` performs the model training for GLMM and GAM to generate the daily resource abundance used in simulation.\
- `02_IBM` contains IBM of optimal foraging fishers under the four different scenario, including `optimal.R`, `time_indifferent.R`, `null_resource.R` and `null.R`.\
- `03_model_perform` handles the calculation of model performance metric and the comparison between different model scenarios.\
\
## Software\
\
```R\
import foobar\
```\
\
## Contributions\
\
Pull requests are welcome. For major changes, please open an issue first\
to discuss what you would like to change.\
\
Please make sure to update tests as appropriate.\
\
## License\
\
[MIT](https://choosealicense.com/licenses/mit/)}