To run the experiments and produce plots as in the paper, you will need python3 and g++ installed, as well as pip, matplotlib, and libboost-dev. You will also need to install the Abseil Python Common Libraries. On Ubuntu Focal for example:

`$ sudo apt install g++ python3 python3-pip python3-matplotlib libboost-dev`

`$ pip3 install absl-py`

To then generate the experimental results and figures in the paper (which a copy has been preserved for you in `experiments`), from the command line type

`$ make figures`

All the experimental results and figures will then live in a new directory called `experiments`.
