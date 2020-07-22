# AC-DC-SBM
Assortative-Constrained Stochastic Block Models

## Run

To run the AC-DC-SBM algorithm in Julia, execute the following command:

`> julia AC_DCSBM.jl 'dataset_name' seed nb_runs`

### Example

`> julia AC_DCSBM.jl 'A-4-05-02-101' 1234 100`

## Data format

**Graph file.** The graph file is a 3-column and m-row file, where m is the number of edges in the graph. Each column is separated by a single space, and each line correponds to one edge in the graph. In a row, the first value is the number of the first sample, the second value is the number of the second sample, and the third value is the weigth of the edge:

|     |     |     |
|-----|-----|-----|
| a<sub>1</sub> | b<sub>1</sub> | w<sub>1</sub> |
| a<sub>2</sub> | b<sub>2</sub> | w<sub>2</sub> |
| ... | ... | ... |
| a<sub>m</sub> | b<sub>m</sub> | w<sub>m</sub> |

**Important**: Graph files must have the `.link` extension. Some graphs are provided within folder `/data` in this repository.

**Label file.** The label file contains the ground-truth community of each sample of the dataset, and is a N-row file, where N is the number of samples. Each line presents the label of the i-th sample:

y<sub>1</sub>

y<sub>2</sub>

...

y<sub>N</sub>

**Important**: Labels files must have the `.label` extension. Some labels are provided within folder `/data` in this repository.
