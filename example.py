from getDistances import *

# Example of running the code

# gap weight: higher = less penalty for gaps, more likely to include alignment with gaps
gap_weight = -1

getCompositionality('ExampleData.csv',"AlignmentWeights.csv", gap_weight, "results.csv")